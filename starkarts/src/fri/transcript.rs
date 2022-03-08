use byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};

use crate::{
    field_elem::FieldElement, Evaluation, FriCommitment, FriOpenings, MerkleError, MerkleTree,
};

pub enum FriError {
    InvalidDomain,
    InvalidIndex,
    NotPowerOfTwoSized,
    InvalidLeaf,
    TooManySamplesRequested,
    InvalidCommitment,
    NotColinear,
}

impl From<MerkleError> for FriError {
    fn from(e: MerkleError) -> Self {
        match e {
            MerkleError::NotPowerOfTwo => Self::NotPowerOfTwoSized,
            MerkleError::IndexTooLarge => Self::InvalidIndex,
            MerkleError::InvalidLeaf => Self::InvalidLeaf,
        }
    }
}
#[derive(Debug)]
pub struct FriProof<F> {
    pub commitment: FriCommitment<F>,
    pub openings: FriOpenings<F>,
}

pub struct Fri<F> {
    pub offset: F,
    pub omega: F,
    pub domain_size: usize,
    // TODO: consider removing
    pub expansion_factor: usize,
    pub num_colinearity_checks: usize,
}
impl<F: FieldElement> Fri<F> {
    /// The number of rounds in the FRI interaction. Note that
    /// the the protocol terminates as soon as the number of colinearity checks
    /// reaches 1/4 of the codeword length. At that point, there are too few possible checks
    /// random subsampling to have any benefit, so we just check the entire codeword.
    ///
    /// TODO: Calculate analytically instead of looping
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

    /// Create a FRI proof that a particular codeword is of low degree
    pub fn prove(&self, codeword: Evaluation<F>) -> Result<FriProof<F>, FriError> {
        if codeword.len() != self.domain_size {
            return Err(FriError::InvalidDomain);
        }

        // Commit to the codewords
        let (commitment, codewords) =
            FriCommitment::commit(codeword, self.num_rounds(), self.omega, self.offset);

        // Get random indices "from the verifier"
        let indices =
            self.determine_sample_indices(&commitment.fiat_shamir(), codewords[0].len())?;

        let openings = FriOpenings::open(&commitment, &indices, &codewords)?;

        Ok(FriProof {
            commitment,
            openings,
        })
    }

    // TODO
    pub fn verify(&self, proof: &FriProof<F>) -> Result<(), FriError> {
        let num_rounds = self.num_rounds();
        let commitment = &proof.commitment;
        let openings = &proof.openings;
        let final_codeword = &commitment.final_codeword;

        // Verify that the final codeword is low degree
        commitment.verify_final_codeword(
            final_codeword.len() / self.expansion_factor,
            self.omega.pow(2 * num_rounds),
        )?;

        // Compute the indices to check for colinearity
        let indices = self.determine_sample_indices(&commitment.fiat_shamir(), self.domain_size)?;

        // Verify that the prover has provided enough openings and roots.
        // Without this check, we might not check every proof because of
        // a Zip terminates when the first iterator runs out.
        if commitment.merkle_roots.len() != openings.0.len() {
            return Err(FriError::InvalidCommitment);
        }
        for (round, ([current_root, next_root], opening)) in commitment
            .merkle_roots
            .array_windows::<2>()
            .zip(openings.0.iter())
            .enumerate()
        {
            let round_indices: Vec<(usize, usize, usize)> = indices
                .iter()
                .map(|&idx| {
                    let reduced_domain_len = self.domain_size >> round + 1;
                    let a_idx = idx % reduced_domain_len;
                    let b_idx = a_idx + reduced_domain_len;
                    (a_idx, b_idx, a_idx)
                })
                .collect();
            let alpha = F::sample(&commitment.partial_fiat_shamir(round));
            opening.verify(
                &round_indices,
                alpha,
                self.omega,
                self.offset,
                current_root,
                next_root,
            )?;
        }

        Ok(())
    }

    #[inline]
    /// Read an index from the random bytes, reducing (mod size)
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

    /// Sample `self.num_colinearity_checks` random indices from the provided seed
    /// using a cryptographic hash function. The returned indices are in the range
    /// 0..length_of_original_codeword/2
    pub fn determine_sample_indices(
        &self,
        seed: &[u8],
        length_of_original_codeword: usize,
    ) -> Result<Vec<usize>, FriError> {
        let mut indices = Vec::with_capacity(self.num_colinearity_checks);
        let mut reduced_indices = Vec::with_capacity(self.num_colinearity_checks);

        let mut sample_seed: Vec<u8> = Vec::with_capacity(seed.len() + 8);
        sample_seed.extend_from_slice(seed);

        let mut iter_number = 0;

        let max_index = length_of_original_codeword / 2;
        let size_after_all_reductions = length_of_original_codeword >> self.num_rounds();

        while indices.len() < self.num_colinearity_checks {
            // Append the sample number to the seed
            sample_seed
                .write_u64::<LittleEndian>(iter_number)
                .expect("write to vec must not fail");

            // Use a hash function to generate a random number, then reduce mod half the size of the initial codeword
            // Recall that colinearity checks in FRI verify that current_codeword[i] and current_codeword[i + len/2] lie on a
            // straight line with next_codeword[i]. Picking indices between 0 and len/2 guarantees that we can safey add len/2
            // to each index without running out of bounds.
            let index = self.sample_index(&MerkleTree::hash(&sample_seed), max_index)?;
            let reduced_index = index % size_after_all_reductions;

            // Verify that the index doesn't collide even when reduced to work for the smallest codeword
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
}
