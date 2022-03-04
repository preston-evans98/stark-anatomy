use blake2::digest::{consts::U32, generic_array::GenericArray};

use crate::{
    field_elem::FieldElement, Evaluation, FriCommitment, FriError, MerkleError, MerkleTree,
};

#[derive(Debug)]
pub struct FriOpening<F>(Vec<Vec<ConsistencyProof<F>>>);

#[derive(Debug)]
pub struct MerkleProof<F> {
    pub merkle_path: Vec<GenericArray<u8, U32>>,
    pub leaf: F,
    pub idx: usize,
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
        Ok(Self {
            merkle_path,
            leaf,
            idx,
        })
    }

    pub fn root(&self) -> Option<&GenericArray<u8, U32>> {
        self.merkle_path.first()
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

    /// Verify the proof by
    pub fn verify(&self, _current_root: &Evaluation<F>, _next_root: &Evaluation<F>) -> bool {
        // check that b.idx = a.idx + current_codeword.len() / 2;

        // Check each merkle proof verifies, using the same root for a and b

        // check that the points are colinear

        todo!()
    }
}

impl<F: FieldElement> FriOpening<F> {
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
    ) -> Result<Vec<ConsistencyProof<F>>, FriError> {
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

        Ok(consistency_proofs)
    }
}
