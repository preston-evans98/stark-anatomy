use crate::{field_elem::FieldElement, GoldilocksField, P as GoldilocksPrime};

/// An element of the Goldilocks field. Not necessarily in canonical form
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Eq, Ord)]
pub struct GoldilocksElement(u64);

impl GoldilocksElement {
    pub fn from_u32(item: u32) -> Self {
        Self(item as u64)
    }

    pub fn from_u64(item: u64) -> Self {
        Self(item)
    }

    pub fn canonical_form(self) -> Self {
        if self.0 >= GoldilocksPrime {
            return Self(self.0 - GoldilocksPrime);
        }
        self
    }

    pub fn to_canonical(&mut self) {
        if self.0 >= GoldilocksPrime {
            self.0 -= GoldilocksPrime
        }
    }
}

impl PartialEq<u64> for GoldilocksElement {
    fn eq(&self, other: &u64) -> bool {
        self.canonical_form().0 == *other
    }
}

impl From<i32> for GoldilocksElement {
    fn from(item: i32) -> Self {
        if item < 0 {
            return Self(GoldilocksPrime - ((-item) as u64));
        }
        Self(item as u64)
    }
}

impl From<u32> for GoldilocksElement {
    fn from(item: u32) -> Self {
        Self(item as u64)
    }
}

impl From<u64> for GoldilocksElement {
    fn from(item: u64) -> Self {
        Self(item)
    }
}

impl std::ops::Add<&Self> for GoldilocksElement {
    type Output = Self;

    fn add(self, rhs: &Self) -> Self::Output {
        Self(GoldilocksField::add(self.0, rhs.0))
    }
}

impl std::ops::Sub<&Self> for GoldilocksElement {
    type Output = Self;

    fn sub(self, rhs: &Self) -> Self::Output {
        Self(GoldilocksField::subtract(self.0, rhs.0))
    }
}

impl std::ops::Mul<&Self> for GoldilocksElement {
    type Output = Self;

    fn mul(self, rhs: &Self) -> Self::Output {
        Self(GoldilocksField::reduce_multiplication(
            self.0.widening_mul(rhs.0),
        ))
    }
}

impl<T: Into<Self>> std::ops::Add<T> for GoldilocksElement {
    type Output = Self;

    fn add(self, rhs: T) -> Self::Output {
        Self(GoldilocksField::add(self.0, rhs.into().0))
    }
}

impl<T: Into<Self>> std::ops::Sub<T> for GoldilocksElement {
    type Output = Self;

    fn sub(self, rhs: T) -> Self::Output {
        Self(GoldilocksField::subtract(self.0, rhs.into().0))
    }
}

impl<T: Into<Self>> std::ops::Mul<T> for GoldilocksElement {
    type Output = Self;

    fn mul(self, rhs: T) -> Self::Output {
        Self(GoldilocksField::reduce_multiplication(
            self.0.widening_mul(rhs.into().0),
        ))
    }
}

impl<T: Into<Self>> std::ops::Div<T> for GoldilocksElement {
    type Output = Self;

    fn div(self, rhs: T) -> Self::Output {
        Self(GoldilocksField::multiply(self.0, rhs.into().inverse().0))
    }
}

impl std::ops::Div<&Self> for GoldilocksElement {
    type Output = Self;

    fn div(self, rhs: &Self) -> Self::Output {
        Self(GoldilocksField::multiply(self.0, rhs.inverse().0))
    }
}

impl std::ops::Neg for GoldilocksElement {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self(GoldilocksField::subtract(GoldilocksPrime, self.0))
    }
}

impl FieldElement for GoldilocksElement {
    fn is_zero(self) -> bool {
        self.0 == 0
    }

    fn zero() -> Self {
        Self(0)
    }

    fn one() -> Self {
        Self(1)
    }

    fn xgcd(self, rhs: Self) -> (Self, Self, Self) {
        let (a, b, c) = GoldilocksField::xgcd(self.0, rhs.0);
        (Self(a), Self(b), Self(c))
    }

    fn inverse(self) -> Self {
        let (a, _, _) = self.xgcd(Self(GoldilocksPrime));
        a
    }

    /// Raises self to the power `exponent`
    fn pow(self, exponent: usize) -> Self {
        pow_iter(self, Self(1), exponent)
    }
}

fn pow_iter(
    value: GoldilocksElement,
    acc: GoldilocksElement,
    exponent: usize,
) -> GoldilocksElement {
    match exponent {
        0 => acc,
        1 => value * acc,
        x if (x & 1 == 0) => pow_iter(value * value, acc, exponent >> 1),
        _ => pow_iter(value.clone(), acc * value, exponent - 1),
    }
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use crate::{field_elem::FieldElement, GoldilocksElement, P};

    #[test]
    fn test_invert() {
        let two: GoldilocksElement = 2.into();
        let result = two.inverse();
        assert_eq!(result, 9223372034707292161)
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(10000))]
        #[test]
        fn adds_correctly(a: u64, b: u64) {
            let g_result = GoldilocksElement::from_u64(a) + GoldilocksElement::from_u64(b);
            let correct_result = ((a as u128) + (b as u128)) % (P as u128);
            prop_assert_eq!((g_result.0 % P) as u128, correct_result)
        }

        #[test]
        fn subtracts_correctly(a: u64, b: u64) {
            let g_result = GoldilocksElement::from_u64(a) - GoldilocksElement::from_u64(b);
            let correct_result = if a > b {
                ((a as u128) - (b as u128)) % (P as u128)
            } else {
                (((a as u128) + (P as u128)) - (b as u128)) % (P as u128)
            };
            prop_assert_eq!((g_result.0 % P) as u128, correct_result)
        }

        #[test]
        fn multiplies_correctly(a: u64, b: u64) {
            let g_result = GoldilocksElement::from_u64(a) * GoldilocksElement::from_u64(b);
            let correct_result =
                ((a as u128) * (b as u128)) % (P as u128)
            ;
            prop_assert_eq!((g_result.0 % P) as u128, correct_result)
        }

        #[test]
        fn inverts_correctly(a: u64) {
            let a = GoldilocksElement::from_u64(a);
            let inverse = a.inverse();
            prop_assert_eq!((a * inverse).0 % P, 1)
        }

        #[test]
        fn divides_correctly(a: u64, b: u64) {
            let div_result = GoldilocksElement::from_u64(a) / GoldilocksElement::from_u64(b);
            let actual = div_result * GoldilocksElement::from_u64(b);
            prop_assert_eq!(actual.0 % P, a)
        }

        #[test]
        fn negates_correctly(a: u64) {
            let a = GoldilocksElement::from_u64(a);
            let neg = -a;
            prop_assert_eq!((a + neg).0 % P , 0)
        }

        #[test]
        fn pows_correctly(a: u64, b: usize) {
            let a = GoldilocksElement::from_u64(a);
            let r  = a.pow(b).0 % P;
            let correct_result = pow_128(a.0 as u128, b, P as u128);
            prop_assert_eq!(r as u128, correct_result)
        }
    }

    fn pow_128(val: u128, exponent: usize, modulus: u128) -> u128 {
        match exponent {
            0 => 1,
            1 => val,
            x if (x % 2) == 0 => pow_128((val * val) % modulus, exponent / 2, modulus),
            x => (val * pow_128(val, exponent - 1, modulus)) % modulus,
        }
    }
}
