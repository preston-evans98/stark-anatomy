use std::ops::{Add, Div, Mul, Neg, Sub};

use crate::uint::Uint;

pub trait FieldElement:
    Copy
    + Add<Self, Output = Self>
    + Sub<Self, Output = Self>
    + Mul<Self, Output = Self>
    + Div<Self, Output = Self>
    + for<'a> Add<&'a Self, Output = Self>
    + for<'a> Sub<&'a Self, Output = Self>
    + for<'a> Mul<&'a Self, Output = Self>
    + for<'a> Div<&'a Self, Output = Self>
    + Neg<Output = Self>
    + PartialOrd
    + PartialEq
    + std::fmt::Debug
    + From<usize>
{
    const byte_length: usize;
    fn xgcd(self, rhs: Self) -> (Self, Self, Self);

    /// Returns the additive identity
    fn zero() -> Self;
    /// Returns the multiplicative identity
    fn one() -> Self;
    fn is_zero(self) -> bool;

    fn inverse(self) -> Self;
    fn pow(self, exponent: usize) -> Self;

    fn zeros(len: usize) -> Vec<Self> {
        std::iter::repeat(Self::zero()).take(len).collect()
    }
    fn to_bytes_le(&self, out: &mut [u8]);
    fn sample(random: &[u8]) -> Self;
}
