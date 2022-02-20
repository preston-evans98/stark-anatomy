use blake2::digest::consts::U32;
use blake2::digest::generic_array::GenericArray;
use blake2::{Blake2b, Digest};

/// A very inefficient implementation of a Merkle tree.
/// This is translated from Python code at https://aszepieniec.github.io/stark-anatomy/basic-tools
/// TODO: Replace this with an efficient (and idiomatic) implementation
pub struct MerkleTree {}
#[derive(Debug)]
pub enum MerkleError {
    NotPowerOfTwo,
    IndexTooLarge,
}

impl MerkleTree {
    pub fn new() -> Self {
        Self {}
    }
    pub fn hash(input: &[u8]) -> GenericArray<u8, U32> {
        let mut hasher = Blake2b::<U32>::new();
        hasher.update(input);
        hasher.finalize()
    }
    fn hash_2(input: [GenericArray<u8, U32>; 2]) -> GenericArray<u8, U32> {
        let mut hasher = Blake2b::<U32>::new();
        hasher.update(input[0]);
        hasher.update(input[1]);
        hasher.finalize()
    }

    fn _commit(leaves: &[GenericArray<u8, U32>]) -> Result<GenericArray<u8, U32>, MerkleError> {
        if !leaves.len().is_power_of_two() {
            return Err(MerkleError::NotPowerOfTwo);
        }
        if leaves.len() == 1 {
            return Ok(leaves[0]);
        }
        let (left, right) = leaves.split_at(leaves.len() / 2);
        Ok(Self::hash_2([Self::_commit(left)?, Self::_commit(right)?]))
    }

    fn _open(
        index: usize,
        leaves: &[GenericArray<u8, U32>],
    ) -> Result<Vec<GenericArray<u8, U32>>, MerkleError> {
        if !leaves.len().is_power_of_two() {
            return Err(MerkleError::NotPowerOfTwo);
        }
        if leaves.len() <= index {
            return Err(MerkleError::IndexTooLarge);
        }
        if leaves.len() == 2 {
            return Ok(vec![leaves[1 - index]]);
        }
        let middle = leaves.len() / 2;
        if index < middle {
            let mut left = Self::_open(index, &leaves[..middle])?;
            let right = Self::_commit(&leaves[middle..])?;
            left.push(right);
            return Ok(left);
        }
        let mut left = Self::_open(index - middle, &leaves[middle..])?;
        let right = Self::_commit(&leaves[..middle])?;
        left.push(right);
        Ok(left)
    }

    fn _verify(
        root: GenericArray<u8, U32>,
        index: usize,
        path: &[GenericArray<u8, U32>],
        leaf: GenericArray<u8, U32>,
    ) -> Result<bool, MerkleError> {
        if index >= path.len() {
            return Err(MerkleError::IndexTooLarge);
        }
        if path.len() == 1 {
            if index == 0 {
                return Ok(root == Self::hash_2([leaf, path[0]]));
            }
            return Ok(root == Self::hash_2([path[0], leaf]));
        } else {
            if index % 2 == 0 {
                return Self::_verify(root, index >> 1, &path[1..], Self::hash_2([leaf, path[0]]));
            }
            Self::_verify(root, index >> 1, &path[1..], Self::hash_2([path[0], leaf]))
        }
    }
    /// Compute the merkle root of an array of data.
    /// The length of the array must be a power of two.
    pub fn commit(leaves: &[&[u8]]) -> Result<GenericArray<u8, U32>, MerkleError> {
        let mut hashed_leaves = Vec::with_capacity(leaves.len());
        for &item in leaves.iter() {
            hashed_leaves.push(Self::hash(item))
        }
        Self::_commit(&hashed_leaves)
    }
    /// Compute the merkle root of an array of hashes of leaves.
    /// The length of the array must be a power of two.
    pub fn commit_preprocessed(
        hashed_leaves: &[GenericArray<u8, U32>],
    ) -> Result<GenericArray<u8, U32>, MerkleError> {
        Self::_commit(hashed_leaves)
    }
    /// Create a merkle proof for the leaf at some index.
    pub fn open(index: usize, leaves: &[&[u8]]) -> Result<Vec<GenericArray<u8, U32>>, MerkleError> {
        let mut hashed_leaves = Vec::with_capacity(leaves.len());
        for &item in leaves.iter() {
            hashed_leaves.push(Self::hash(item))
        }
        Self::_open(index, &hashed_leaves)
    }
    /// Verify a merkle proof
    pub fn verify(
        root: GenericArray<u8, U32>,
        index: usize,
        path: &[GenericArray<u8, U32>],
        leaf: &[u8],
    ) -> Result<bool, MerkleError> {
        let leaf = Self::hash(leaf);
        Self::_verify(root, index, path, leaf)
    }
}

#[test]
fn test_multiple_update() {
    let input = [[1u8; 32], [1u8; 32]];
    let mut hasher = Blake2b::<U32>::new();
    hasher.update(input[0]);
    hasher.update(input[1]);
    let output = hasher.finalize();

    let mut hasher2 = Blake2b::<U32>::new();
    hasher2.update(&[1u8; 64]);
    let output2 = hasher2.finalize();
    assert_eq!(output, output2)
}
