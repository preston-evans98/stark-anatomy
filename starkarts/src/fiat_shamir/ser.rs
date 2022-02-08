/// A canonical serialization of objects. This is important for security.
/// If the serialization of each object is not unique, an attacker can
/// manipulate the randomness output by fiat-shamir and break soundness.
///
/// Note: the tutorial at https://aszepieniec.github.io/stark-anatomy/basic-tools
/// uses Python's `Pickle` library, which does not guarantee this property.
/// Pickle may provide a small enough set of possible serializations for each object
/// that this soundess cannot be broken - further Googling is needed to confirm this.
pub trait CanonicalSer {
    fn canon_serialize_to_vec(&self) -> Vec<u8> {
        let mut out = Vec::new();
        self.canon_serialize(&mut out);
        out
    }
    fn canon_serialize(&self, _out: &mut Vec<u8>) {}
}
