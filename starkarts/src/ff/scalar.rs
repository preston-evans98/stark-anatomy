use std::ops::{Add, Div, Mul, Rem, Sub, Neg, AddAssign};

use num::{BigUint,  Zero, One, Integer};

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

    fn pow(&self, exponent: usize) -> Self {
        pow_iter(self.clone(), Self::one(), exponent)
    }

    
}

fn  pow_iter<S: Scalar>(value: S, acc: S, exponent: usize) -> S {
    match exponent {
        0 => S::one(),
        1 => value * acc,
        x if x.is_even() => {
            pow_iter(value.clone() * value, acc, exponent >> 1)
        }
        _  => pow_iter(value.clone(), acc*value, exponent - 1),
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

#[cfg(test)]
mod tests {
    use num::{BigUint, Num};

    use crate::{scalar::Scalar, DefaultField};


    #[test]
    fn pow_test() {
        let elem = BigUint::from(12345u16);
        let result = elem.pow(2345) %
        BigUint::from_str_radix("270497897142230380135924736767050121217", 10).unwrap();
        assert_eq!(
            result,
            BigUint::from_str_radix("23520667101221308275390517182423299422", 10).unwrap()
        )
    }
}
