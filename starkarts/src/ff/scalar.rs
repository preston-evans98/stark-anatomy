use std::ops::{Add, Div, Mul, Rem, Sub, Neg};

use num::{BigInt,  Zero, One, Integer};

pub trait SimpleNum: 
    Clone + 
    std::cmp::PartialOrd<Self> 
    + Sized + Mul<Output = Self>
    + Sub<Output = Self>
    + Add<Output = Self>
    + Rem<Output = Self>
    + std::fmt::Debug {
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
}

/// All instances of any field implement the trait Scalar. 
/// The definition of a field is as follows:
/// 1. It has an addition operation that is associative and commutative, plus an additive identity.
/// 2. Every element of the field has an additive inverse, and the field is closed under addition
/// 3. It has a multiplication operation that is associative and commutative, plus a multiplicative identity
/// 4. Every element of the field has a multiplicative inverse, and the field is closed under multiplication
/// 5. Multiplication distributes over addition: (a + b) * c = a*c + b*c
pub trait Scalar: SimpleNum
    // Copy
	
    + Div<Output = Self>
    
{
    
    /// Implements the extended euclidean algorithm for finding greatest common denominators
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

    /// Returns the multiplicative inverse of element
    fn inverse(&self) -> Self {
        Self::one() / self.clone()
    }
}

fn  pow_iter<S: SimpleNum>(value: S, acc: S, exponent: usize) -> S {
    match exponent {
        0 => S::one(),
        1 => value * acc,
        x if x.is_even() => {
            pow_iter(value.clone() * value, acc, exponent >> 1)
        }
        _  => pow_iter(value.clone(), acc*value, exponent - 1),
    }
}

impl SimpleNum for BigInt {
	fn is_zero(&self) -> bool {
		*self == <Self as SimpleNum>::zero()
	}

	fn zero() -> Self {
		<BigInt as Zero>::zero()
	}

	fn one() -> Self {
		<BigInt as One>::one()
	}
}

impl Scalar for BigInt{}

pub trait CopyScalar: Scalar + Copy{}

pub trait SignedScalar: CopyScalar + Neg<Output = Self>{}

#[cfg(test)]
mod tests {
    use num::{BigInt, Num, FromPrimitive, Integer};

    use crate::{scalar::Scalar, DefaultField};

    #[test]
    fn gcd_test() {
        let a = BigInt::from_str_radix("270497897142230380135924736767050121217", 10).unwrap(); 
        let b = BigInt::from_i128(123128739127).unwrap();

        assert_eq!(b.extended_gcd(&a).gcd, BigInt::from_i128(1).unwrap());
        // assert_eq!(a.extended_gcd(&b).gcd, BigInt::from_i128(1).unwrap());
    }

    #[test]
    fn gcd_test_2() {
        let a = BigInt::from_str_radix("270497897142230380135924736767050121217", 10).unwrap(); 
        let b = BigInt::from_str_radix("270497897142230380135924736767050121216", 10).unwrap(); 

        // assert_eq!(b.extended_gcd(&a).gcd, BigInt::from_i128(1).unwrap());
        assert_eq!(Scalar::extended_gcd(b.clone(), a.clone()).0, b);
        // assert_eq!(a.extended_gcd(&b).gcd, BigInt::from_i128(1).unwrap());
    }

    #[test]
    fn pow_test() {
        let elem = BigInt::from(12345u16);
        let result = elem.pow(2345) %
        BigInt::from_str_radix("270497897142230380135924736767050121217", 10).unwrap();
        assert_eq!(
            result,
            BigInt::from_str_radix("23520667101221308275390517182423299422", 10).unwrap()
        )
    }
}
