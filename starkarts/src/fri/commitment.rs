use blake2::digest::{
    consts::U32, core_api::CoreWrapper, generic_array::GenericArray, ExtendableOutput, Update,
    XofReader,
};
use sha3::{Shake256, Shake256Core};

use crate::{field_elem::FieldElement, CanonicalSer, Evaluation, MerkleTree};

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

    pub fn prover_fiat_shamir(&self) -> [u8; 32] {
        let mut hasher = Shake256::default();
        self.fiat_shamir_update(&mut hasher);
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
