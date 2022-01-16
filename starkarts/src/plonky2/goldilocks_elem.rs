use crate::{
    scalar::{CopyScalar, Scalar, SimpleNum},
    GoldilocksField, GoldilocksMath, P as GoldilocksPrime,
};

/// An element of the Goldilocks field. Not necessarily in canonical form
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Eq, Ord)]
pub struct GoldilocksElement(u64);

impl GoldilocksElement {
    pub fn from_u32(item: u32) -> Self {
        Self(item as u64)
    }
}

impl From<u32> for GoldilocksElement {
    fn from(item: u32) -> Self {
        Self(item as u64)
    }
}

impl std::ops::Add<Self> for GoldilocksElement {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self(GoldilocksField::add(self.0, rhs.0))
    }
}

impl std::ops::Sub<Self> for GoldilocksElement {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self(GoldilocksField::subtract(self.0, rhs.0))
    }
}

impl std::ops::Mul<Self> for GoldilocksElement {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self(GoldilocksField::reduce_multiplication(
            self.0.widening_mul(rhs.0),
        ))
    }
}

impl std::ops::Div<Self> for GoldilocksElement {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        Self(GoldilocksField::multiply(self.0, rhs.inverse().0))
    }
}

impl std::ops::Rem<Self> for GoldilocksElement {
    type Output = Self;

    fn rem(self, rhs: Self) -> Self::Output {
        // let (gcd, x, y) = GoldilocksField::xgcd(self.0, rhs.0)
        unimplemented!()
    }
}

impl SimpleNum for GoldilocksElement {
    fn is_zero(&self) -> bool {
        self.0 == 0
    }

    fn zero() -> Self {
        Self(0)
    }

    fn one() -> Self {
        Self(1)
    }

    fn zeros(len: usize) -> Vec<Self> {
        std::iter::repeat(Self::zero()).take(len).collect()
    }

    fn pow(&self, exponent: usize) -> Self {
        pow_iter(*self, Self::one(), exponent as u64)
    }
}
impl Scalar for GoldilocksElement {
    fn extended_gcd(x: Self, y: Self) -> (Self, Self, Self) {
        let (a, b, c) = GoldilocksField::xgcd(x.0, y.0);
        (Self(a), Self(b), Self(c))
    }

    fn inverse(&self) -> Self {
        let (a, _, _) = self.xgcd(Self(GoldilocksPrime));
        a
    }
}

impl CopyScalar for GoldilocksElement {}

impl GoldilocksMath for GoldilocksElement {
    fn xgcd(self, rhs: Self) -> (Self, Self, Self) {
        let (a, b, c) = GoldilocksField::xgcd(self.0, rhs.0);
        (Self(a), Self(b), Self(c))
    }

    fn inverse(self) -> Self {
        let (a, _, _) = self.xgcd(Self(GoldilocksPrime));
        a
    }

    /// Raises self to the power `exponent`
    fn pow(self, exponent: u64) -> Self {
        pow_iter(self, Self(1), exponent)
    }

    // fn zeros(len: usize) -> Vec<Self> {
    //     std::iter::repeat(Self(0)).take(len).collect()
    // }
    // fn zero() -> Self {
    //     Self(0)
    // }

    // fn one() -> Self {
    //     Self(1)
    // }
}

fn pow_iter(value: GoldilocksElement, acc: GoldilocksElement, exponent: u64) -> GoldilocksElement {
    match exponent {
        0 => GoldilocksElement(0),
        1 => value * acc,
        x if (x & 1 == 0) => pow_iter(value * value, acc, exponent >> 1),
        _ => pow_iter(value.clone(), acc * value, exponent - 1),
    }
}
