use crate::{
    scalar::{Scalar, SignedScalar},
    DivByZeroError,
};

#[derive(Debug, Clone, PartialEq)]
/// A polynomial, represented as a list of coefficients, with the lowest degree coefficient first.
/// For example, 5x^2 + 1 is represented [1, 0, 5]
pub struct Polynomial<S> {
    coefficients: Vec<S>,
}

impl<S> Polynomial<S>
where
    S: Scalar,
{
    pub fn from(coefs: Vec<S>) -> Self {
        let mut p = Self {
            coefficients: coefs,
        };
        p._trim_degree();
        p
    }

    pub fn from_ref(coefs: &[S]) -> Self {
        let mut p = Self {
            coefficients: coefs.to_vec(),
        };
        p._trim_degree();
        p
    }
    /// Returns the degree of the polynomial.
    /// The degree of the zero polynomial is defined to be -1
    pub fn degree(&self) -> isize {
        self.coefficients
            .len()
            .try_into()
            .expect("Must be able to convert degree of any constructed polynomial to isize")
    }

    pub fn zero() -> Self {
        Self {
            coefficients: Vec::new(),
        }
    }

    pub fn is_zero(&self) -> bool {
        self.degree() == -1
    }

    pub fn leading_coefficient(&self) -> Option<&S> {
        self.coefficients.last()
    }

    fn _trim_degree(&mut self) {
        while self.coefficients.ends_with(&[S::zero()]) {
            self.coefficients.pop();
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
        if self.is_zero() {
            return rhs;
        } else if rhs.is_zero() {
            return self;
        }

        let (long_iter, short_iter) = if self.coefficients.len() > rhs.coefficients.len() {
            (self.coefficients.into_iter(), rhs.coefficients.into_iter())
        } else {
            (rhs.coefficients.into_iter(), self.coefficients.into_iter())
        };

        return Self::from(
            long_iter
                .zip(short_iter.chain(std::iter::repeat(S::zero())))
                .map(|(l, r)| l + r)
                .collect(),
        );
    }
}

impl<S> std::ops::Sub for Polynomial<S>
where
    S: Scalar,
{
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        if self.is_zero() {
            return rhs;
        } else if rhs.is_zero() {
            return self;
        }

        let (long_iter, short_iter) = if self.coefficients.len() > rhs.coefficients.len() {
            (self.coefficients.into_iter(), rhs.coefficients.into_iter())
        } else {
            (rhs.coefficients.into_iter(), self.coefficients.into_iter())
        };

        return Self::from(
            long_iter
                .zip(short_iter.chain(std::iter::repeat(S::zero())))
                .map(|(l, r)| l - r)
                .collect(),
        );
    }
}

impl<S> std::ops::Mul for Polynomial<S>
where
    S: Scalar,
{
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        if self.is_zero() || rhs.is_zero() {
            return Self::zero();
        }
        let mut coefs = S::zeros((self.degree() + rhs.degree() + 1) as usize);

        for (i, c1) in self.coefficients.iter().enumerate() {
            // If c1 is zero then c1 * c2 = 0 * c2 = 0, so we can skip this iteration
            if c1.is_zero() {
                continue;
            }
            for (j, c2) in rhs.coefficients.iter().enumerate() {
                coefs[i + j] += c1.clone() * c2.clone()
            }
        }

        return Polynomial::from(coefs);
    }
}

// TODO: for better performance, we can use FFTs instead of
// polynomial long divison. For now, this will do.
impl<S> std::ops::Div for Polynomial<S>
where
    S: Scalar,
{
    type Output = Result<(Self, Self), DivByZeroError>;

    fn div(self, divisor: Self) -> Self::Output {
        if divisor.is_zero() {
            return Err(DivByZeroError);
        }

        if self.degree() < divisor.degree() {
            return Ok((Self::zero(), self));
        }

        let result_degree_upper_bound = (self.degree() - divisor.degree() + 1) as usize;

        let mut remainder = self;
        let mut quotient_coefs = S::zeros(result_degree_upper_bound);
        for _ in 0..result_degree_upper_bound {
            if remainder.degree() < divisor.degree() {
                break;
            }
            let coefficient = remainder
                .leading_coefficient()
                .expect("we just checked that the remainder is non-zero")
                .clone()
                - divisor
                    .leading_coefficient()
                    .expect("divisor must be non-zero")
                    .clone();

            let shift = (remainder.degree() - divisor.degree()) as usize;
            let subtrahend = Polynomial::from(
                std::iter::repeat(S::zero())
                    .take(shift)
                    .chain(std::iter::once(coefficient.clone()))
                    .collect(),
            ) * divisor.clone();
            quotient_coefs[shift] = coefficient;
            remainder = remainder - subtrahend;
        }
        let quotient = Polynomial::from(quotient_coefs);

        // TODO: remove
        return Ok((quotient, remainder));
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
        let a = Polynomial::from_ref(&[1, 2, 3, 4, 5]);
        let b = Polynomial::from_ref(&[1, 2, 3]);

        assert_eq!(
            (a + b).coefficients,
            Polynomial::from_ref(&[2, 4, 6, 4, 5]).coefficients
        )
    }
}
