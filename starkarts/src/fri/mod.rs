use blake2::digest::{consts::U32, generic_array::GenericArray};
use byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};

use crate::{
    field_elem::FieldElement, CanonObjectTag, CanonicalSer, Evaluation, MerkleError, MerkleTree,
    ProofStream,
};

#[derive(Debug, Copy, Clone)]
pub enum FriError {
    InvalidDomain,
    /// Trying to sample too many random indices for too small a codeword
    TooManySamplesRequested,
    NotPowerOfTwoSized,
    InvalidIndex,
}

impl From<MerkleError> for FriError {
    fn from(e: MerkleError) -> Self {
        match e {
            MerkleError::NotPowerOfTwo => Self::NotPowerOfTwoSized,
            MerkleError::IndexTooLarge => Self::InvalidIndex,
        }
    }
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
        codeword: Evaluation<F>,
        proof_stream: &mut ProofStream,
    ) -> Result<Vec<usize>, FriError> {
        if codeword.len() != self.domain_size {
            return Err(FriError::InvalidDomain);
        }

        // Commit to the codewords
        let codewords = self.commit(codeword, proof_stream);
        // Get random indices "from the verifier"
        let indices = self.sample_indices(
            &proof_stream.prover_fiat_shamir(),
            codewords[1].len(),
            codewords[codewords.len() - 1].len(),
            self.num_colinearity_checks,
        )?;

        // Open the commitment at verifier's indices
        let mut indices_this_round = indices.clone();
        for [cur, next] in codewords.array_windows::<2>() {
            indices_this_round
                .iter_mut()
                .for_each(|idx| *idx = *idx % cur.len() / 2);

            self.query(cur, next, &indices_this_round, proof_stream)?;
        }
        return Ok(indices);
    }

    pub fn commit(
        &self,
        mut codeword: Evaluation<F>,
        proof_stream: &mut ProofStream,
    ) -> Vec<Evaluation<F>> {
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
        proof_stream.push(Box::new(codeword.clone()));
        codewords.push(codeword);

        codewords
    }

    #[inline]
    /// Read an index from the random bytes, reducing (mod size)/
    ///
    /// In this implementation, we assume that size is always a power of two
    /// so there is no need to read the high-order-bytes.
    pub fn sample_index(&self, bytes: &[u8], size: usize) -> Result<usize, FriError> {
        if !size.is_power_of_two() {
            return Err(FriError::NotPowerOfTwoSized);
        }
        let mut cursor = std::io::Cursor::new(bytes);
        let raw_index = cursor
            .read_u64::<LittleEndian>()
            .expect("can always read 8 bytes from a 32 byte array");
        return Ok((raw_index % (size as u64)) as usize);
    }

    pub fn sample_indices(
        &self,
        seed: &[u8],
        size: usize,
        reduced_size: usize,
        num_samples: usize,
    ) -> Result<Vec<usize>, FriError> {
        if num_samples >= 2 * reduced_size {
            return Err(FriError::TooManySamplesRequested);
        }

        let mut indices = Vec::with_capacity(num_samples);
        let mut reduced_indices = Vec::with_capacity(num_samples);

        let mut sample_seed: Vec<u8> = Vec::with_capacity(seed.len() + 8);
        sample_seed.extend_from_slice(seed);
        let mut iter_number = 0;
        while indices.len() < num_samples {
            // Append the sample number to the seed
            sample_seed
                .write_u64::<LittleEndian>(iter_number)
                .expect("write to vec must not fail");

            // Use a hash function to generate a random number
            let index = self.sample_index(&MerkleTree::hash(&sample_seed), size)?;
            let reduced_index = index % reduced_size;

            // Verify that the number hasn't collided with another index
            if !reduced_indices.contains(&reduced_index) {
                reduced_indices.push(reduced_index);
                indices.push(index);
            }

            // Truncate the sample_number from the seed
            sample_seed.truncate(seed.len());
            iter_number += 1;
        }

        Ok(indices)
    }

    pub fn query(
        &self,
        current_codeword: &Evaluation<F>,
        next_codeword: &Evaluation<F>,
        c_indices: &[usize],
        proof_stream: &mut ProofStream,
    ) -> Result<Vec<usize>, FriError> {
        let a_indices = c_indices.clone();
        let b_indices: Vec<usize> = c_indices
            .iter()
            .map(|&idx| idx + (current_codeword.len() / 2))
            .collect();

        for ((&a, &b), &c) in a_indices.iter().zip(b_indices.iter()).zip(c_indices.iter()) {
            proof_stream.push(Box::new((
                current_codeword.at(a),
                current_codeword.at(b),
                next_codeword.at(c),
            )))
        }

        let current_preprocessed = current_codeword.merkle_preprocess();
        let next_preprocessed = &next_codeword.merkle_preprocess();

        for ((&a, &b), &c) in a_indices.iter().zip(b_indices.iter()).zip(c_indices.iter()) {
            proof_stream.push(Box::new(MerkleTree::open_preprocessed(
                a,
                &current_preprocessed,
            )?));
            proof_stream.push(Box::new(MerkleTree::open_preprocessed(
                b,
                &current_preprocessed,
            )?));
            proof_stream.push(Box::new(MerkleTree::open_preprocessed(
                c,
                &next_preprocessed,
            )?));
        }

        let mut result = a_indices.to_vec();
        result.extend(b_indices.into_iter());

        Ok(result)
    }
}

impl CanonicalSer for GenericArray<u8, U32> {
    fn canon_serialize(&self, out: &mut Vec<u8>) {
        out.push(CanonObjectTag::Bytes32 as u8);
        out.extend(self[..].iter())
    }
}

impl<S: CanonicalSer> CanonicalSer for Vec<S> {
    fn canon_serialize(&self, out: &mut Vec<u8>) {
        out.push(CanonObjectTag::Vec as u8);
        out.write_u64::<LittleEndian>(self.len() as u64)
            .expect("serialization to vec must succeed");
        for item in self.iter() {
            item.canon_serialize(out)
        }
    }
}

impl<F: CanonicalSer> CanonicalSer for (F, F, F) {
    fn canon_serialize(&self, out: &mut Vec<u8>) {
        out.push(CanonObjectTag::Triple as u8);
        self.0.canon_serialize(out);
        self.1.canon_serialize(out);
        self.2.canon_serialize(out);
    }
}
