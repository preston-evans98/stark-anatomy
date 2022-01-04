use num::BigUint;
use starkarts::{DefaultField, PrimeField};

pub fn main() {
    dbg!(DefaultField::p());
    dbg!(BigUint::from_bytes_be(&[1, 2, 3, 4, 5, 6, 7, 8, 1]));
}
