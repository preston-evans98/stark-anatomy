use crate::CanonicalSer;
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
}

// impl CanonicalSer for ProofStream {
//     fn canon_serialize(&self, mut out: &mut Vec<u8>) {
//         for item in self.objects.iter() {
//             item.canon_serialize(&mut out);
//         }
//     }
// }

impl CanonicalSer for ProofStream {
    fn canon_serialize(&self, mut out: &mut Vec<u8>) {
        for item in self.objects.iter() {
            item.canon_serialize(&mut out);
        }
    }
}

use std::rc::Rc;

use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake128,
};

// use hex_literal::hex;

// let mut hasher = Shake128::default();
// hasher.update(b"abc");
// let mut reader = hasher.finalize_xof();
// let mut res1 = [0u8; 10];
// reader.read(&mut res1);
// assert_eq!(res1, hex!("5881092dd818bf5cf8a3"));
