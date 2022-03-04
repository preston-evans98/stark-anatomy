use crate::{field_elem::FieldElement, CanonObjectTag, CanonicalSer, Evaluation};
use byteorder::{LittleEndian, WriteBytesExt};
use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake256,
};

#[derive(Debug)]
/// A ProofStream contains all of the messages that would have
/// been sent between a Prover and Verifier if the proof were
/// interactive.
pub struct ProofStream<F, S: CanonicalSer> {
    pub objects: Vec<S>,
    pub codeword: Option<Evaluation<F>>,
    pub cursor: usize,
}
impl<F: FieldElement, S: CanonicalSer> ProofStream<F, S> {
    pub fn new() -> Self {
        Self {
            objects: vec![],
            cursor: 0,
            codeword: None,
        }
    }
    pub fn push(&mut self, obj: S) {
        self.objects.push(obj)
    }
    pub fn advance_cursor(&mut self) {
        self.cursor += 1
    }
    pub fn num_items(&self) -> usize {
        self.objects.len() + if let Some(_) = self.codeword { 1 } else { 0 }
    }

    // pub fn next<'a, 's>(&'s mut self) -> Option<&'a Box<dyn CanonicalSer>>
    // where
    //     's: 'a,
    // {
    //     if self.cursor < self.objects.len() {
    //         self.cursor += 1;
    //         return Some(&self.objects[self.cursor - 1]);
    //     }
    //     None
    // }
    // pub fn next<'a, 's>(&'s mut self, out: &mut Vec<&'a Box<dyn CanonicalSer>>) -> Option<()>
    // where
    //     's: 'a,
    // {
    //     if self.cursor < self.objects.len() {
    //         self.cursor += 1;
    //         out.push(&self.objects[self.cursor - 1]);
    //         return Some(());
    //     }
    //     None
    // }

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

        if self.cursor > self.objects.len() {
            if let Some(ref codeword) = self.codeword {
                codeword.canon_serialize(&mut out);
            }
        }
    }
}

impl<F: FieldElement, S: CanonicalSer> CanonicalSer for ProofStream<F, S> {
    fn canon_serialize(&self, mut out: &mut Vec<u8>) {
        out.write_u8(CanonObjectTag::ProofStream as u8)
            .expect("write to vec must succeed");
        out.write_u32::<LittleEndian>(self.num_items() as u32)
            .expect("write to vec must succeed");
        for item in self.objects.iter() {
            item.canon_serialize(&mut out);
        }
    }
}
