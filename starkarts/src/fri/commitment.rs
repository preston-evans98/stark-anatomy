use blake2::digest::{
    consts::U32, core_api::CoreWrapper, generic_array::GenericArray, ExtendableOutput, Update,
    XofReader,
};
use sha3::{Shake256, Shake256Core};

use crate::{field_elem::FieldElement, CanonicalSer, Evaluation, FriError, MerkleTree, Polynomial};

#[derive(Debug)]
/// An efficiently queryable commitment to a low degree polynomial.
pub struct FriCommitment<F> {
    /// The merkle root of each successive codeword, ordered from largest to smallest.
    pub merkle_roots: Vec<GenericArray<u8, U32>>,
    /// The final codeword, whose size should be smaller than the original by a factor of 2^num_rounds
    pub final_codeword: Evaluation<F>,
}

impl<F: FieldElement> FriCommitment<F> {
    /// Use FRI to commit to a codeword (remember, a codeword is just a polynomial in evaluation form).
    /// Use the finite cyclic group with generator `omega` as the evaluation domain, but start from
    /// `omega * offset` instead of `omega`. Use `num_rounds` rounds of split-and-fold to reduce the size of the codeword
    /// to a manageable size. The final codeword will have shrunk by a factor of 2^(num_rounds).
    ///
    /// Returns a FriProof and Vec all of the intermediate codewords created while generating the proof
    pub fn commit(
        mut codeword: Evaluation<F>,
        num_rounds: usize,
        mut omega: F,
        mut offset: F,
    ) -> (Self, Vec<Evaluation<F>>) {
        let mut merkle_roots = Vec::with_capacity(num_rounds + 1);
        let mut intermediate_codewords = Vec::with_capacity(num_rounds);
        for r in 0..num_rounds {
            // Commit to the polynomial
            let root = MerkleTree::commit_preprocessed(&codeword.merkle_preprocess())
                .expect("Preprocessed input must have power-of-two length");
            merkle_roots.push(root);

            // Prepare for the next round. If there are no rounds left, exit
            if r == num_rounds - 1 {
                break;
            }

            // Use the fiat-shamir heuristic to get a "random" alpha
            // from the "verifier"
            let alpha = F::sample(&prover_fiat_shamir(&merkle_roots));

            // Halve the degree of the codeword
            let new_codeword = codeword.split_and_fold(alpha, omega, offset);

            // Collect the old codeword
            intermediate_codewords.push(codeword);
            codeword = new_codeword;

            omega = omega * omega;
            offset = offset * offset;
        }

        // Add the last codeword to the list of intermediate results
        intermediate_codewords.push(codeword.clone());

        // Create the FRI proof
        let proof = Self {
            merkle_roots,
            final_codeword: codeword,
        };

        // Return the proof and intermediate computations
        (proof, intermediate_codewords)
    }

    /// Compute the hash of the entire commitment
    pub fn fiat_shamir(&self) -> [u8; 32] {
        let mut hasher = Shake256::default();
        self.fiat_shamir_update(&mut hasher);
        let mut reader = hasher.finalize_xof();
        let mut result = [0u8; 32];
        reader.read(&mut result);
        result
    }

    /// Compute the hash of the commitment up to (but not including) some index
    pub fn partial_fiat_shamir(&self, up_to: usize) -> [u8; 32] {
        let mut hasher = Shake256::default();
        for root in self.merkle_roots.iter().take(up_to) {
            hasher.update(&root[..]);
        }
        let mut reader = hasher.finalize_xof();
        let mut result = [0u8; 32];
        reader.read(&mut result);
        result
    }

    /// Update the state of `hasher` with the (entire) contents of this commitment
    pub fn fiat_shamir_update(&self, hasher: &mut CoreWrapper<Shake256Core>) {
        for root in self.merkle_roots.iter() {
            hasher.update(&root[..]);
        }
        hasher.update(&self.final_codeword.canon_serialize_to_vec());
    }

    /// Verifies that the `final_codeword` is low degree and matches the merkle commitment
    /// given the generator `omega` of the cyclic group which should have been used to compute the *final* codeword
    /// and the expected degree
    pub fn verify_final_codeword(
        &self,
        expected_degree: usize,
        final_omega: F,
    ) -> Result<(), FriError> {
        // Verify that the final_codeword matches the last merkle root
        match self.merkle_roots.last() {
            Some(root) => {
                if MerkleTree::commit_preprocessed(&self.final_codeword.merkle_preprocess())?[..]
                    != root[..]
                {
                    return Err(FriError::InvalidCommitment);
                }
            }
            None => {
                return Err(FriError::InvalidCommitment);
            }
        }

        // Verify that omega generates a cyclic subgroup with the correct order
        let degree = self.final_codeword.len() + 1;
        if final_omega.inverse() != final_omega.pow(degree) {
            return Err(FriError::InvalidDomain);
        }

        // Compute the domain consisting of [omega^i for i in range 0..degree+1]
        let mut domain = Vec::with_capacity(degree + 1);
        let mut omega = final_omega;
        for _ in 0..degree + 1 {
            domain.push(omega);
            omega = omega * final_omega;
        }

        // Interpolate a polynomial from the evaluations and verify its degree
        let interpolated =
            Polynomial::<F>::interpolate_domain(&domain, &self.final_codeword.evaluations)
                .map_err(|_| FriError::InvalidDomain)?;
        if interpolated.degree() > expected_degree as isize {
            return Err(FriError::InvalidCommitment);
        }

        Ok(())
    }
}

/// Run the fiat shamir heuristic over the proof stream so far
fn prover_fiat_shamir(input: &Vec<GenericArray<u8, U32>>) -> [u8; 32] {
    let mut hasher = Shake256::default();
    for root in input.iter() {
        hasher.update(&root[..]);
    }
    let mut reader = hasher.finalize_xof();
    let mut result = [0u8; 32];
    reader.read(&mut result);
    result
}
