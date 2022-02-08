use crate::{CanonObjectTag, CanonicalSer};
use byteorder::{LittleEndian, WriteBytesExt};
use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake256,
};

/// A ProofStream contains all of the messages that would have
/// been sent between a Prover and Verifier if the proof were
/// interactive.
pub struct ProofStream {
    objects: Vec<Box<dyn CanonicalSer>>,
    cursor: usize,
}
impl ProofStream {
    pub fn new() -> Self {
        Self {
            objects: vec![],
            cursor: 0,
        }
    }
    pub fn push(&mut self, obj: Box<dyn CanonicalSer>) {
        self.objects.push(obj)
    }
    pub fn next(&mut self) -> Option<&Box<dyn CanonicalSer>> {
        if self.cursor < self.objects.len() {
            self.cursor += 1;
            return Some(&self.objects[self.cursor - 1]);
        }
        None
    }

    pub fn prover_fiat_shamir(&self) -> [u8; 32] {
        let mut hasher = Shake256::default();
        hasher.update(&self.canon_serialize_to_vec());
        let mut reader = hasher.finalize_xof();
        let mut result = [0u8; 32];
        reader.read(&mut result);
        result
    }

    pub fn verifier_fiat_shamir(&self) -> [u8; 32] {
        let mut serialized = Vec::new();
        let mut hasher = Shake256::default();
        self.serialize_up_to_cursor(&mut serialized);
        hasher.update(&mut serialized);

        let mut reader = hasher.finalize_xof();
        let mut result = [0u8; 32];
        reader.read(&mut result);
        result
    }

    fn serialize_up_to_cursor(&self, mut out: &mut Vec<u8>) {
        out.write_u8(CanonObjectTag::ProofStream as u8)
            .expect("write to vec must succeed");
        out.write_u32::<LittleEndian>(self.cursor as u32)
            .expect("write to vec must succeed");
        for item in self.objects.iter().take(self.cursor) {
            item.canon_serialize(&mut out);
        }
    }
}

impl CanonicalSer for ProofStream {
    fn canon_serialize(&self, mut out: &mut Vec<u8>) {
        out.write_u8(CanonObjectTag::ProofStream as u8)
            .expect("write to vec must succeed");
        out.write_u32::<LittleEndian>(self.objects.len() as u32)
            .expect("write to vec must succeed");
        for item in self.objects.iter() {
            item.canon_serialize(&mut out);
        }
    }
}
