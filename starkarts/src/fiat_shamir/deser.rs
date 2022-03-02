use byteorder::{LittleEndian, ReadBytesExt};
use core::slice::Iter;

use crate::{field_elem::FieldElement, CanonicalSer, ProofStream};

// use crate::CanonSerializableObject;

impl CanonObjectTag {
    pub fn deser<R: std::io::Read, F: FieldElement + 'static>(
        cursor: &mut R,
    ) -> Result<Box<dyn CanonicalSer>, DeserializationError> {
        let tag = cursor.read_u8()?;
        Ok(Box::new(match tag.try_into()? {
            CanonObjectTag::ProofStream => Self::read_proof_stream::<R, F>(cursor)?,
            CanonObjectTag::Bytes32 => todo!(),
            CanonObjectTag::Evaluation => todo!(),
            CanonObjectTag::Triple => todo!(),
            CanonObjectTag::Vec => todo!(),
        }))
    }

    fn read_proof_stream<R: std::io::Read, F: FieldElement + 'static>(
        cursor: &mut R,
    ) -> Result<ProofStream<F>, DeserializationError> {
        let len = cursor.read_u32::<LittleEndian>()?;
        let mut proof = ProofStream::new();
        for _ in 0..len {
            let item = Self::deser::<R, F>(cursor)?;
            proof.push(item);
        }
        Ok(proof)
    }
}

pub trait CanonDeser: Sized {
    fn canon_deser(iter: Iter<u8>) -> Result<Self, DeserializationError>;
}

pub enum DeserializationError {
    MissingObjectTag,
    MalformedData,
}

impl From<std::io::Error> for DeserializationError {
    fn from(_e: std::io::Error) -> Self {
        Self::MalformedData
    }
}

#[repr(u8)]
pub enum CanonObjectTag {
    ProofStream = 0,
    Bytes32 = 1,
    Evaluation = 2,
    Triple = 3,
    Vec = 4,
}

impl TryFrom<u8> for CanonObjectTag {
    type Error = DeserializationError;

    fn try_from(value: u8) -> Result<Self, Self::Error> {
        match value {
            0 => Ok(Self::ProofStream),
            1 => Ok(Self::Bytes32),
            2 => Ok(Self::Evaluation),
            3 => Ok(Self::Triple),
            4 => Ok(Self::Vec),
            _ => Err(DeserializationError::MissingObjectTag),
        }
    }
}
