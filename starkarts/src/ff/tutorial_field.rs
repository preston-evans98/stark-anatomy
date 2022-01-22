use bigint::U256;

use crate::{
    field::{NoNthRootError, PrimeField},
    uint::Uint,
};

// use std::ops::{Add, Div, Mul, Rem, Sub};

// use crate::scalar::SimpleNum;

// use super::scalar::Scalar;

lazy_static::lazy_static! {
    static ref P: U256 = U256::from_dec_str("270497897142230380135924736767050121217").unwrap();
    static ref GENERATOR: U256 = U256::from_dec_str("85408008396924667383611388730472331217").unwrap();
    static ref ORDER: U256 = U256::from_dec_str("664613997892457936451903530140172288").unwrap();
}

#[derive(Debug, Clone)]
pub struct StarkTutorialField;
impl PrimeField for StarkTutorialField {
    type Elem = U256;
    #[inline]
    fn p() -> U256 {
        *P
    }

    #[inline]
    fn generator() -> U256 {
        *GENERATOR
    }

    fn primitive_nth_root(n: U256) -> Result<U256, NoNthRootError> {
        if n >= *ORDER || n.is_odd() {
            return Err(NoNthRootError);
        }
        let mut root = ORDER.clone();
        let mut order = ORDER.clone();
        while order != n {
            root = root * root;
            order = order >> 1;
        }
        Ok(root)
    }

    fn sample(random_bytes: &[u8]) -> U256 {
        U256::from_big_endian(random_bytes) % Self::p()
    }

    fn reduce(n: Self::Elem) -> Self::Elem {
        n % Self::p()
    }
}

#[cfg(test)]
mod tests {
    use bigint::U256;

    use crate::field::PrimeField;

    use super::StarkTutorialField;

    #[test]
    fn test_subtract() {}

    #[test]
    fn test_inverse() {
        let one: U256 = 1.into();
        let a: U256 = StarkTutorialField::p() - one;
        let expected_result =
            U256::from_dec_str("270497897142230380135924736767050121216").unwrap();
        assert_eq!(a, expected_result);

        assert_eq!(StarkTutorialField::inverse(a), expected_result);
    }
}
