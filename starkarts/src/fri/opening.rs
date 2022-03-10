use blake2::digest::{consts::U32, generic_array::GenericArray};

use crate::{
    field_elem::FieldElement, CanonicalSer, Evaluation, FriCommitment, FriError, MerkleError,
    MerkleTree, Polynomial,
};

#[derive(Debug)]
/// All of the witnesses needed to verify a Fri commitment.
pub struct FriOpenings<F>(pub Vec<FriRoundOpening<F>>);

#[derive(Debug)]
/// The witnesses needed to verify a single round of FRI
pub struct FriRoundOpening<F>(pub Vec<ConsistencyProof<F>>);

impl<F: FieldElement> FriRoundOpening<F> {
    // TODO: double check
    /// Verifies a single round of FRI. In other words, verifies that two merkle roots `current` and `next`
    /// both commit to polynomials, and that `next` is related to `current` by the FRI equation.

    pub fn verify(
        &self,
        indices: &[(usize, usize, usize)],
        alpha: F,
        omega: F,
        offset: F,
        current_root: &GenericArray<u8, U32>,
        next_root: &GenericArray<u8, U32>,
    ) -> Result<(), FriError> {
        // Be sure to use every set of indices. If there aren't enough
        // Consistency proofs to do so, the Commitment was invalid.
        if self.0.len() != indices.len() {
            return Err(FriError::InvalidCommitment);
        }
        for (proof, indices) in self.0.iter().zip(indices.iter()) {
            let a_x = offset * omega.pow(indices.0);
            let b_x = offset * omega.pow(indices.1);
            proof.verify(a_x, b_x, alpha, indices, current_root, next_root)?;
        }
        Ok(())
    }
}

#[derive(Debug)]
pub struct MerkleProof<F> {
    pub merkle_path: Vec<GenericArray<u8, U32>>,
    pub leaf: F,
}

impl<F: FieldElement> MerkleProof<F> {
    pub fn new(
        leaf: F,
        idx: usize,
        hashed_leaves: &[GenericArray<u8, U32>],
    ) -> Result<Self, MerkleError> {
        if MerkleTree::hash(&leaf.canon_serialize_to_vec()) != hashed_leaves[idx] {
            return Err(MerkleError::InvalidLeaf);
        }
        let merkle_path = MerkleTree::open_preprocessed(idx, hashed_leaves)?;
        Ok(Self { merkle_path, leaf })
    }

    pub fn root(&self) -> Option<&GenericArray<u8, U32>> {
        self.merkle_path.first()
    }

    pub fn verify(&self, root: &GenericArray<u8, U32>, idx: usize) -> Result<(), MerkleError> {
        match MerkleTree::verify(
            root.clone(),
            idx,
            &self.merkle_path,
            &self.leaf.canon_serialize_to_vec(),
        ) {
            Ok(true) => Ok(()),
            _ => Err(MerkleError::InvalidLeaf),
        }
    }
}

#[derive(Debug)]
/// A proof that two Reed Solomon codewords are equivalent at two indices.
/// To verify, check that 1.
/// the merkle proofs hold, that A and C are part of the s
pub struct ConsistencyProof<F> {
    pub a: MerkleProof<F>,
    pub b: MerkleProof<F>,
    pub c: MerkleProof<F>,
}

impl<F: FieldElement> ConsistencyProof<F> {
    pub fn new(a: MerkleProof<F>, b: MerkleProof<F>, c: MerkleProof<F>) -> Self {
        Self { a, b, c }
    }

    /// Verify the proof by checking the merkle roots and verifying that the
    /// three points are colinear.
    pub fn verify(
        &self,
        a_x: F,
        b_x: F,
        alpha: F,
        indices: &(usize, usize, usize),
        current_root: &GenericArray<u8, U32>,
        next_root: &GenericArray<u8, U32>,
    ) -> Result<(), FriError> {
        let (a_idx, b_idx, c_idx) = indices;
        // Check the merkle proofs
        self.a.verify(current_root, *a_idx)?;
        self.b.verify(current_root, *b_idx)?;
        self.c.verify(next_root, *c_idx)?;

        // check that the points are colinear
        if !Polynomial::are_colinear(vec![
            (a_x, self.a.leaf),
            (b_x, self.b.leaf),
            (alpha, self.c.leaf),
        ]) {
            return Err(FriError::NotColinear);
        }
        Ok(())
    }
}

impl<F: FieldElement> FriOpenings<F> {
    pub fn open(
        commitment: &FriCommitment<F>,
        at: &[usize],
        codewords: &Vec<Evaluation<F>>,
    ) -> Result<Self, FriError> {
        let num_rounds = commitment.merkle_roots.len();
        let mut openings = Vec::with_capacity(num_rounds);
        let mut indices = vec![];
        indices.extend_from_slice(at);
        for [cur, next] in codewords.array_windows::<2>() {
            indices
                .iter_mut()
                .for_each(|idx| *idx = *idx % (cur.len() / 2));

            openings.push(Self::query(cur, next, &indices)?)
        }

        Ok(Self(openings))
    }

    pub fn query(
        current_codeword: &Evaluation<F>,
        next_codeword: &Evaluation<F>,
        indices: &[usize],
    ) -> Result<FriRoundOpening<F>, FriError> {
        let b_offset = current_codeword.len() / 2;
        let mut consistency_proofs = Vec::with_capacity(indices.len());

        // TODO: this is horribly inefficient. Reimplement merkle tree with intermediate results
        let current_preprocessed = current_codeword.merkle_preprocess();
        let next_preprocessed = next_codeword.merkle_preprocess();

        for &idx in indices.iter() {
            let b_idx = idx + b_offset;
            consistency_proofs.push(ConsistencyProof::new(
                MerkleProof::new(current_codeword.at(idx), idx, &current_preprocessed)?,
                MerkleProof::new(current_codeword.at(b_idx), b_idx, &current_preprocessed)?,
                MerkleProof::new(next_codeword.at(idx), idx, &next_preprocessed)?,
            ))
        }

        Ok(FriRoundOpening(consistency_proofs))
    }
}
