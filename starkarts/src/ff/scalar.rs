use std::ops::{Add, Div, Mul, Rem, Sub, Neg, AddAssign};

use num::{BigUint,  Zero, One};

pub trait Scalar:
    // Copy
	Clone + 
    std::cmp::PartialOrd<Self> 
	+ Sized
    + Div<Output = Self>
    + Mul<Output = Self>
    + Sub<Output = Self>
    + Add<Output = Self>
    + AddAssign
    + Rem<Output = Self>
    + std::fmt::Debug
{
    fn is_zero(&self) -> bool;
    fn zero() -> Self;
    fn one() -> Self;
    fn extended_gcd(x: Self, y: Self) -> (Self, Self, Self) {
        let (mut old_r, mut r) = (x, y);
        let (mut old_s, mut s) = (Self::one(), Self::zero());
        let (mut old_t, mut t) = (Self::zero(), Self::one());

        while !r.is_zero() {
            let quotient = old_r.clone() / r.clone();
            (old_r, r) = (r.clone(), old_r - quotient.clone() * r);
            (old_s, s) = (s.clone(), old_s - quotient.clone() * s);
            (old_t, t) = (t.clone(), old_t - quotient * t);
        }
        (old_s, old_t, old_r)
    }

    fn zeros(len: usize) -> Vec<Self> {
        std::iter::repeat(Self::zero()).take(len).collect()
    }
}

impl Scalar for BigUint {
	fn is_zero(&self) -> bool {
		*self == <Self as Scalar>::zero()
	}

	fn zero() -> Self {
		<BigUint as Zero>::zero()
	}

	fn one() -> Self {
		<BigUint as One>::one()
	}
}

pub trait SignedScalar: Scalar + Neg<Output = Self>{}
