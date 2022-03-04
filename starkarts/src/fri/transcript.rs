use byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};

use crate::{
    field_elem::FieldElement, Evaluation, FriCommitment, FriOpening, MerkleError, MerkleTree,
};

pub enum FriError {
    InvalidDomain,
    InvalidIndex,
    NotPowerOfTwoSized,
    InvalidLeaf,
    TooManySamplesRequested,
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
    pub opening: FriOpening<F>,
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
        let indices = self.sample_indices(
            &commitment.prover_fiat_shamir(),
            codewords[1].len(),
            codewords.last().expect("codewords must not be empty").len(),
            self.num_colinearity_checks,
        )?;

        let opening = FriOpening::open(&commitment, &indices, &codewords)?;

        Ok(FriProof {
            commitment,
            opening,
        })
    }

    // TODO
    pub fn verify() {}

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

    /// Sample `num_samples` random indices from the provided seed
    ///  using a cryptographic hash function
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
}
