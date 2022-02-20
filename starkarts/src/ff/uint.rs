use std::ops::{Add, Div, Mul, Sub};

use bigint::U256;

pub trait Uint:
    Copy
    + Add<Self, Output = Self>
    + Sub<Self, Output = Self>
    + Mul<Self, Output = Self>
    + Div<Self, Output = Self>
    + PartialOrd
    + PartialEq
    + std::fmt::Debug
{
    /// Returns true if the elment is zero (the additive identity)
    fn is_zero(&self) -> bool;
    /// Returns the additive identity (0)
    fn zero() -> Self;
    /// Returns the multiplicative identity (1)
    fn one() -> Self;
    fn zeros(len: usize) -> Vec<Self> {
        std::iter::repeat(Self::zero()).take(len).collect()
    }

    /// Raises self to the power `exponent`
    fn pow(&self, exponent: usize) -> Self {
        pow_iter(self.clone(), Self::one(), exponent)
    }

    fn is_even(&self) -> bool;
    fn is_odd(&self) -> bool;
    fn to_bytes_le(&self, out: &mut [u8]);
}

/// All instances of any field implement the trait Scalar.
/// The definition of a field is as follows:
/// 1. It has an addition operation that is associative and commutative, plus an additive identity.
/// 2. Every element of the field has an additive inverse, and the field is closed under addition
/// 3. It has a multiplication operation that is associative and commutative, plus a multiplicative identity
/// 4. Every element of the field has a multiplicative inverse, and the field is closed under multiplication
/// 5. Multiplication distributes over addition: (a + b) * c = a*c + b*c

fn pow_iter<S: Uint>(value: S, acc: S, exponent: usize) -> S {
    match exponent {
        0 => S::one(),
        1 => value * acc,
        x if (x & 1) == 0 => pow_iter(value.clone() * value, acc, exponent >> 1),
        _ => pow_iter(value.clone(), acc * value, exponent - 1),
    }
}

impl Uint for U256 {
    fn is_zero(&self) -> bool {
        self.is_zero()
    }

    fn zero() -> Self {
        <U256>::zero()
    }

    fn one() -> Self {
        <U256>::one()
    }

    fn is_even(&self) -> bool {
        !self.is_odd()
    }
    fn is_odd(&self) -> bool {
        self.bit(0)
    }
    fn to_bytes_le(&self, out: &mut [u8]) {
        self.to_little_endian(out);
    }
}
