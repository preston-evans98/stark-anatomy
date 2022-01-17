use std::ops::{Add, Div, Mul, Neg, Sub};

pub trait FieldElement:
    Copy
    + Add<Self, Output = Self>
    + Sub<Self, Output = Self>
    + Mul<Self, Output = Self>
    + Div<Self, Output = Self>
    + Neg<Output = Self>
    + PartialOrd
    + PartialEq
    + std::fmt::Debug
{
    fn xgcd(self, rhs: Self) -> (Self, Self, Self);

    // Returns the additive identity
    fn zero() -> Self;
    // Returns the multiplicative identity
    fn one() -> Self;
    fn is_zero(self) -> bool;

    fn inverse(self) -> Self;
    fn pow(self, exponent: usize) -> Self;

    fn zeros(len: usize) -> Vec<Self> {
        std::iter::repeat(Self::zero()).take(len).collect()
    }
}
