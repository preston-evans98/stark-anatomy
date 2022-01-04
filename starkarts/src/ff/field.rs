use num::{BigUint, Integer, Num};

// use std::ops::{Add, Div, Mul, Rem, Sub};

use super::field_elem::FieldElement;

pub struct DivByZeroError;
pub struct NoNthRootError;

lazy_static::lazy_static! {
    static ref P: BigUint = BigUint::from_str_radix("270497897142230380135924736767050121217", 10).unwrap();
    static ref GENERATOR: BigUint = BigUint::from_str_radix("85408008396924667383611388730472331217", 10).unwrap();
    static ref ORDER: BigUint = BigUint::from_str_radix("664613997892457936451903530140172288", 10).unwrap();
}

pub trait PrimeField {
    type Elem: FieldElement;
    fn zero() -> Self::Elem {
        <Self::Elem as FieldElement>::zero()
    }
    fn one() -> Self::Elem {
        <Self::Elem as FieldElement>::one()
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

pub struct DefaultField;
impl PrimeField for DefaultField {
    type Elem = BigUint;
    fn p() -> Self::Elem {
        P.clone()
    }

    fn generator() -> Self::Elem {
        GENERATOR.clone()
    }

    fn primitive_nth_root(n: Self::Elem) -> Result<Self::Elem, NoNthRootError> {
        if n >= ORDER.clone() || n.is_odd() {
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
        BigUint::from_bytes_be(random_bytes) % Self::p()
    }
}
