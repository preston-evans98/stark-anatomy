use std::ops::{Add, Div, Mul, Rem, Sub};

pub struct DivByZeroError;
pub trait FieldElement:
    Copy
    + std::cmp::PartialOrd<Self>
    + Div<Output = Self>
    + Mul<Output = Self>
    + Sub<Output = Self>
    + Add<Output = Self>
    + Rem<Output = Self>
{
    fn is_zero(&self) -> bool;
    fn zero() -> Self;
    fn one() -> Self;
    fn extended_gcd(x: Self, y: Self) -> (Self, Self, Self) {
        let (mut old_r, mut r) = (x, y);
        let (mut old_s, mut s) = (Self::one(), Self::zero());
        let (mut old_t, mut t) = (Self::zero(), Self::one());

        while !r.is_zero() {
            let quotient = old_r / r;
            (old_r, r) = (r, old_r - quotient * r);
            (old_s, s) = (s, old_s - quotient * s);
            (old_t, t) = (t, old_t - quotient * t);
        }
        (old_s, old_t, old_r)
    }
}

pub trait PrimeField<E: FieldElement> {
    fn zero() -> E {
        E::zero()
    }
    fn one() -> E {
        E::one()
    }

    fn p() -> E;

    fn add(l: E, r: E) -> E {
        (l + r) % Self::p()
    }

    fn subtract(l: E, r: E) -> E {
        ((Self::p() + l) - r) % Self::p()
    }

    fn multiply(l: E, r: E) -> E {
        (l * r) % Self::p()
    }

    fn negate(value: E) -> E {
        Self::p() - value
    }

    fn inverse(value: E) -> E {
        let (a, _, _) = E::extended_gcd(value, Self::p());
        a
    }

    fn divide(l: E, r: E) -> Result<E, DivByZeroError> {
        if r.is_zero() {
            return Err(DivByZeroError);
        }
        let (a, _, _) = E::extended_gcd(r, Self::p());
        Ok((l * a) % Self::p())
    }
}
