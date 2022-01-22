use crate::uint::Uint;

#[derive(Debug)]
pub struct DivByZeroError;
#[derive(Debug)]
pub struct NoNthRootError;

/// A Finite Field implemented modulo a prime.
pub trait PrimeField: std::fmt::Debug + Clone {
    type Elem: Uint;

    /// Reduce an element to its such form
    fn reduce(n: Self::Elem) -> Self::Elem;

    fn p() -> Self::Elem;
    fn generator() -> Self::Elem;
    fn primitive_nth_root(n: Self::Elem) -> Result<Self::Elem, NoNthRootError>;
    fn sample(random: &[u8]) -> Self::Elem;

    fn zero() -> Self::Elem {
        <Self::Elem as Uint>::zero()
    }
    fn one() -> Self::Elem {
        <Self::Elem as Uint>::one()
    }

    /// Add two field elements. The default implementation assumes that...
    /// 1. Both elements are in the range [0, P] and
    /// 2. Addition of such elements cannot overflow.
    ///
    /// This default can be overridden to relax those assumptions.
    fn add(l: Self::Elem, r: Self::Elem) -> Self::Elem {
        Self::reduce(l + r)
    }

    /// Subtract two field elements. The default implementation assumes that...
    /// 1. Both elements are in the range [0, P] and
    /// 2. Addition of such elements cannot overflow
    ///
    /// This default can be overridden to relax those assumptions.
    fn subtract(l: Self::Elem, r: Self::Elem) -> Self::Elem {
        Self::reduce((Self::p() + l) - r)
    }

    /// Multiply two field elements. The default implementation assumes that...
    /// 1. Both elements are in the range [0, P] and
    /// 2. multiplication of such elements cannot overflow
    ///
    /// This default can be overridden to relax those assumptions.
    fn multiply(l: Self::Elem, r: Self::Elem) -> Self::Elem {
        Self::reduce(l * r)
    }

    /// Negate a field element. The default implementation assumes that...
    /// 1. The element is in the range [0, P]
    ///
    /// This default can be overridden to relax that assumptions.
    fn negate(value: Self::Elem) -> Self::Elem {
        Self::p() - value
    }

    /// Use the extended euclidean algorithm to compute the Greatest Common Denominators on two elements.
    /// The default implementation assumes that...
    /// 1. The elements are in the range [0, P]
    /// 2. Operations on such elements cannot overflow
    /// 3. Operations on canonical elements return canonical elements.
    ///
    /// This default can be overridden to relax that assumptions.
    fn xgcd(lhs: Self::Elem, rhs: Self::Elem) -> (Self::Elem, Self::Elem, Self::Elem) {
        let (mut old_r, mut r) = (lhs, rhs);
        let (mut old_s, mut s) = (Self::one(), Self::zero());
        let (mut old_t, mut t) = (Self::zero(), Self::one());

        while !r.is_zero() {
            let quotient = old_r / r;
            (old_r, r) = (r, Self::subtract(old_r, Self::multiply(quotient, r)));
            (old_s, s) = (s, Self::subtract(old_s, Self::multiply(quotient, s)));
            (old_t, t) = (t, Self::subtract(old_t, Self::multiply(quotient, t)));
        }
        (old_s, old_t, old_r)
    }

    /// Compute the multiplicative inverse of value
    fn inverse(value: Self::Elem) -> Self::Elem {
        let (a, _, _) = Self::xgcd(value, Self::p());
        a
    }

    /// Divide two field elements.
    fn divide(l: Self::Elem, r: Self::Elem) -> Result<Self::Elem, DivByZeroError> {
        if r.is_zero() {
            return Err(DivByZeroError);
        }

        Ok(Self::reduce(l * Self::inverse(r)))
    }

    fn pow(mut base: Self::Elem, mut exponent: usize) -> Self::Elem {
        let mut acc = Self::Elem::one();
        while exponent != 0 {
            if exponent & 1 == 0 {
                base = <Self as PrimeField>::multiply(base, base);
                exponent >>= 1
            } else {
                acc = <Self as PrimeField>::multiply(base, acc);
                exponent -= 1
            }
        }
        acc
    }
}
