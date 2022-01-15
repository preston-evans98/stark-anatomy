use num::{BigInt, Integer, Num};

// use std::ops::{Add, Div, Mul, Rem, Sub};

use crate::scalar::SimpleNum;

use super::scalar::Scalar;

#[derive(Debug)]
pub struct DivByZeroError;
#[derive(Debug)]
pub struct NoNthRootError;

lazy_static::lazy_static! {
    static ref P: BigInt = BigInt::from_str_radix("270497897142230380135924736767050121217", 10).unwrap();
    static ref GENERATOR: BigInt = BigInt::from_str_radix("85408008396924667383611388730472331217", 10).unwrap();
    static ref ORDER: BigInt = BigInt::from_str_radix("664613997892457936451903530140172288", 10).unwrap();
}

pub trait PrimeField: std::fmt::Debug + Clone {
    type Elem: Scalar;
    fn zero() -> Self::Elem {
        <Self::Elem as SimpleNum>::zero()
    }
    fn one() -> Self::Elem {
        <Self::Elem as SimpleNum>::one()
    }

    fn p() -> Self::Elem;
    fn generator() -> Self::Elem;
    fn primitive_nth_root(n: Self::Elem) -> Result<Self::Elem, NoNthRootError>;
    fn sample(random: &[u8]) -> Self::Elem;

    fn add(l: Self::Elem, r: Self::Elem) -> Self::Elem {
        (l + r) % Self::p()
    }

    fn subtract(l: Self::Elem, r: Self::Elem) -> Self::Elem {
        ((Self::p() + l) - r) % Self::p()
    }

    fn multiply(l: Self::Elem, r: Self::Elem) -> Self::Elem {
        (l * r) % Self::p()
    }

    fn negate(value: Self::Elem) -> Self::Elem {
        Self::p() - value
    }

    /// Computes the multiplicative inverse of value
    fn inverse(value: Self::Elem) -> Self::Elem {
        let (a, _, _) = Self::Elem::extended_gcd(value, Self::p());
        a
    }

    fn divide(l: Self::Elem, r: Self::Elem) -> Result<Self::Elem, DivByZeroError> {
        if r.is_zero() {
            return Err(DivByZeroError);
        }
        let (a, _, _) = Self::Elem::extended_gcd(r, Self::p());
        Ok((l * a) % Self::p())
    }
}

#[derive(Debug, Clone)]
pub struct DefaultField;
impl PrimeField for DefaultField {
    type Elem = BigInt;
    fn p() -> Self::Elem {
        P.clone()
    }

    fn generator() -> Self::Elem {
        GENERATOR.clone()
    }

    fn primitive_nth_root(n: Self::Elem) -> Result<Self::Elem, NoNthRootError> {
        if n >= *ORDER || n.is_odd() {
            return Err(NoNthRootError);
        }
        let mut root = ORDER.clone();
        let mut order = ORDER.clone();
        while order != n {
            root = root.pow(2);
            order /= 2usize;
        }
        Ok(root)
    }

    fn sample(random_bytes: &[u8]) -> Self::Elem {
        BigInt::from_bytes_be(num::bigint::Sign::Plus, random_bytes) % Self::p()
    }
}

#[cfg(test)]
mod tests {
    use num::{BigInt, Num};

    use crate::PrimeField;

    use super::DefaultField;

    #[test]
    fn test_inverse() {
        let one: BigInt = 1.into();
        let a: BigInt = DefaultField::p() - one;
        assert_eq!(
            a.clone(),
            BigInt::from_str_radix("270497897142230380135924736767050121216", 10).unwrap()
        );
        assert_eq!(DefaultField::inverse(a), 1.into());
    }
}
