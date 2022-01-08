use crate::scalar::{Scalar, SignedScalar};

#[derive(Debug, Clone)]
/// A polynomial, represented as a list of coefficients.
/// The coefficient of the highest degree (non-zero) term is at index 0.
/// Note that the list of coefficients may have trailing zeros.
pub struct Polynomial<S> {
    coefficients: Vec<S>,
}

impl<S> Polynomial<S>
where
    S: Scalar,
{
    pub fn from(coefs: &[S]) -> Self {
        Self {
            coefficients: coefs.to_vec(),
        }
    }
    /// Returns the degree of the polynomial. Since polynomials may have trailing zeros
    /// in our implementation, this is an O(degree) operation. The degree of the zero
    /// polynomial is defined to be -1
    pub fn degree(&self) -> isize {
        match self
            .coefficients
            .iter()
            .enumerate()
            .filter(|(_idx, coef)| !coef.is_zero())
            .last()
        {
            Some((idx, _coef)) => idx
                .try_into()
                .expect("Must be able to convert degree of any constructed polynomial to isize"),
            None => -1,
        }
    }
}

impl<S> std::ops::Neg for Polynomial<S>
where
    S: SignedScalar,
{
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self {
            coefficients: self.coefficients.iter().map(|c| -c.clone()).collect(),
        }
    }
}

impl<S> std::ops::Add for Polynomial<S>
where
    S: Scalar,
{
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        if self.degree() == -1 {
            return rhs;
        } else if rhs.degree() == -1 {
            return self;
        }
        let (long_iter, short_iter) = if self.coefficients.len() > rhs.coefficients.len() {
            (self.coefficients.iter(), rhs.coefficients.iter())
        } else {
            (rhs.coefficients.iter(), self.coefficients.iter())
        };

        return Self {
            coefficients: long_iter
                .zip(short_iter.chain(std::iter::repeat(&S::zero())))
                .map(|(l, r)| l.clone() + r.clone())
                .collect(),
        };
    }
}

#[cfg(test)]
mod tests {
    use crate::scalar::Scalar;

    use super::Polynomial;

    impl Scalar for i32 {
        fn is_zero(&self) -> bool {
            *self == 0
        }

        fn zero() -> Self {
            0
        }

        fn one() -> Self {
            1
        }
    }

    #[test]
    fn test_add() {
        let a = Polynomial::from(&[1, 2, 3, 4, 5]);
        let b = Polynomial::from(&[1, 2, 3]);

        assert_eq!(
            (a + b).coefficients,
            Polynomial::from(&[2, 4, 6, 4, 5]).coefficients
        )
    }
}
