use blake2::digest::{consts::U32, generic_array::GenericArray};
use byteorder::{LittleEndian, WriteBytesExt};

use crate::{field_elem::FieldElement, CanonObjectTag, CanonicalSer, MerkleTree};

/// A polynomial represented in evaluation form (as a list of points)
/// [f(0), f(1), ..., f(n)]
///
/// Also known as a Reed-Solomon codeword
pub struct Evaluation<F> {
    evaluations: Vec<F>,
}

impl<F: FieldElement> Evaluation<F> {
    /// Creates a related polynomial with half the degree of this one
    /// using the split_and_fold technique from FRI.
    ///
    /// Takes three parameters
    /// 1. alpha: a random number
    /// 1. omega: a generator of a cyclic group with N elements where N is the number of points in the codeword
    /// 1. offset: a random number between 1 and N
    pub fn split_and_fold(&self, alpha: F, omega: F, offset: F) -> Self {
        let mut output = Self {
            evaluations: Vec::with_capacity((self.evaluations.len() + 1) / 2),
        };
        let one = F::one();
        let two: F = 2.into();
        let one_half = two.inverse();
        let mut omega_acc = omega * offset * one_half;
        let midpoint = self.evaluations.len() / 2;
        // Because omega is a cyclic group with order N, omega ^ N = 1
        // Thus (omega ^ (N^2)) ^ 2 = 1, so omega ^ (N^2) is the square root of 1.
        // But since omega generates a cyclic group of order N, omega ^ (N^2) cannot equal 1. If it
        // did, the group would only have order N/2.
        // Therefore omega ^ (N^2) is the square root of one which is not equal to 1. In other words
        // omega ^ (N^2) is -1.
        //
        // We use this fact to easily generate the evaluations f(i) and f(-i)
        // f(i) is just the entry at index `i`, and f(-i) is the entry at index (i + N/2)
        let positive_iter = self.evaluations[..midpoint].into_iter();
        let negative_iter = self.evaluations[midpoint..].into_iter();

        for (&pos, &neg) in positive_iter.zip(negative_iter) {
            let a_over_o = alpha / omega_acc;
            let left = pos * (one + a_over_o);
            let right = neg * (one - a_over_o);
            output.evaluations.push(left + right);
            omega_acc = omega_acc * omega;
        }

        output
    }

    // Returns a vector of the hashes of each coefficient of the polynomial.
    // The vector is padded with hashes of zero, and so is guaranteed to have length 2^N for some N
    pub fn merkle_preprocess<'a>(&'a self) -> Vec<GenericArray<u8, U32>> {
        let mut output = Vec::with_capacity(self.evaluations.len().next_power_of_two());
        // Precompute the hash of zero
        let mut bytes = Vec::with_capacity(F::byte_length);
        F::zero().to_bytes_le(&mut bytes);
        let zero_hash = MerkleTree::hash(&bytes);
        for coef in self.evaluations.iter() {
            coef.to_bytes_le(&mut bytes);
            output.push(MerkleTree::hash(&bytes));
        }
        // hash(Zero) extend the vector to have the correct length
        for _ in self.evaluations.len()..output.capacity() {
            output.push(zero_hash.clone())
        }
        output
    }
}

impl<F: FieldElement> CanonicalSer for Evaluation<F> {
    fn canon_serialize_to_vec(&self) -> Vec<u8> {
        let mut out = Vec::new();
        self.canon_serialize(&mut out);
        out
    }

    fn canon_serialize(&self, out: &mut Vec<u8>) {
        out.push(CanonObjectTag::Evaluation as u8);
        out.write_u64::<LittleEndian>(self.evaluations.len() as u64)
            .expect("serialization to vec cannot fail");
        let mut serialized_eval = Vec::with_capacity(F::byte_length);
        for e in self.evaluations.iter() {
            e.to_bytes_le(&mut serialized_eval);
            out.extend_from_slice(&serialized_eval);
        }
    }
}
