use blake2::digest::{consts::U32, generic_array::GenericArray};

use crate::{
    field_elem::FieldElement, CanonObjectTag, CanonicalSer, Evaluation, MerkleError, MerkleTree,
    Polynomial, ProofStream,
};

pub enum FriError {
    InvalidDomain,
}
pub struct Fri<F> {
    pub offset: F,
    pub omega: F,
    pub domain_size: usize,
    pub expansion_factor: usize,
    pub num_colinearity_checks: usize,
}
impl<F: FieldElement + 'static> Fri<F> {
    /// The number of rounds in the FRI interaction. Note that
    /// the the protocol terminates as soon as the number of colinearity checks
    /// reaches 1/4 of the codeword length. At that point, there are too few possible checks
    /// random subsampling to have any benefit, so we just check the entire codeword.
    ///
    /// TODO: Calculate directly instead of looping
    pub fn num_rounds(&self) -> usize {
        let mut codeword_length = self.domain_size;
        let mut num_rounds = 0;
        while codeword_length > self.expansion_factor
            && 4 * self.num_colinearity_checks < codeword_length
        {
            codeword_length = codeword_length / 2;
            num_rounds += 1;
        }
        num_rounds
    }

    /// Return the evaluation domain
    pub fn eval_domain(&self) -> Vec<F> {
        let mut domain = Vec::with_capacity(self.domain_size);
        for i in 0..self.domain_size {
            // TODO: Make O(n) instead of O(n log n) by using omega as an accumulator
            domain.push(self.offset * self.omega.pow(i));
        }
        domain
    }

    pub fn prove(
        &self,
        codeword: Polynomial<F>,
        _proof_stream: ProofStream,
    ) -> Result<(), FriError> {
        if (codeword.degree() + 1) as usize != self.domain_size {
            return Err(FriError::InvalidDomain);
        }

        todo!()
    }

    pub fn commit(
        &self,
        mut codeword: Evaluation<F>,
        proof_stream: &mut ProofStream,
        round_index: usize,
    ) {
        let mut omega = self.omega;
        let mut offset = self.offset;
        let mut codewords = Vec::with_capacity(self.num_rounds());
        for r in 0..self.num_rounds() {
            // Commit to the polynomial
            let root = MerkleTree::commit_preprocessed(&codeword.merkle_preprocess())
                .expect("Preprocessed input must have power-of-two length");
            proof_stream.push(Box::new(root));

            // Prepare for the next round. If there are no rounds left, exit
            if r == self.num_rounds() - 1 {
                break;
            }

            // Use the fiat-shamir heuristic to get a "random" alpha
            // from the "verifier"
            let alpha = F::sample(&proof_stream.prover_fiat_shamir());

            // Halve the degree of the codeword
            let new_codeword = codeword.split_and_fold(alpha, omega, offset);

            // Collect the old codeword
            codewords.push(codeword);
            codeword = new_codeword;

            omega = omega * omega;
            offset = offset * offset;
        }

        // Add the last codeword to the proof stream
        proof_stream.push(Box::new(codeword));
    }
}

impl CanonicalSer for GenericArray<u8, U32> {
    // fn canon_serialize_to_vec(&self) -> Vec<u8> {
    //     let mut out = Vec::new();
    //     self.canon_serialize(&mut out);
    //     out
    // }

    fn canon_serialize(&self, out: &mut Vec<u8>) {
        out.push(CanonObjectTag::Bytes32 as u8);
        out.extend(self[..].iter())
    }
}
