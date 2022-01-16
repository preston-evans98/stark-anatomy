use std::cmp::Ordering;

use crate::{
    scalar::{CopyScalar, SignedScalar, SimpleNum},
    DivByZeroError, GoldilocksMath,
};

#[derive(Debug)]
pub enum InterpolationError {
    MismatchedCoordinateLists,
    NoPointsProvided,
}

#[derive(Debug, Clone, PartialEq)]
/// A polynomial, represented as a list of coefficients, with the lowest degree coefficient first.
/// For example, 5x^2 + 1 is represented [1, 0, 5]
pub struct Polynomial<S> {
    coefficients: Vec<S>,
}

impl<S: CopyScalar> SimpleNum for Polynomial<S> {
    fn is_zero(&self) -> bool {
        self.degree() == -1
    }

    fn zero() -> Self {
        Self {
            coefficients: vec![],
        }
    }

    fn one() -> Self {
        Self {
            coefficients: vec![S::one()],
        }
    }
}

impl<S> Polynomial<S>
where
    S: CopyScalar,
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
        let len: isize = self
            .coefficients
            .len()
            .try_into()
            .expect("Must be able to convert degree of any constructed polynomial to isize");
        len - 1
    }

    /// Returns the coefficient of the highest degree term
    pub fn leading_coefficient(&self) -> Option<&S> {
        self.coefficients.last()
    }

    /// Uses polynomial long divison to divide `self` by `divisor`.
    /// Returns the quotient and remainder.
    pub fn long_divide(self, divisor: Self) -> Result<(Self, Self), DivByZeroError> {
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
            let coefficient = *remainder
                .leading_coefficient()
                .expect("we just checked that the remainder is non-zero")
                - *(divisor
                    .leading_coefficient()
                    .expect("divisor must be non-zero"));

            let shift = (remainder.degree() - divisor.degree()) as usize;
            let subtrahend = Polynomial::from(
                std::iter::repeat(S::zero())
                    .take(shift)
                    .chain(std::iter::once(coefficient))
                    .collect(),
            ) * &divisor;
            quotient_coefs[shift] = coefficient;
            remainder = remainder - subtrahend;
        }
        let quotient = Polynomial::from(quotient_coefs);
        return Ok((quotient, remainder));
    }

    /// A naive implementation of evaluation at a point
    pub fn evaluate(&self, point: S) -> S {
        println!("Starting evaluation");
        let mut x_i = S::one();
        let mut value = S::zero();
        println!("Looping over {} elements", self.degree() + 1);
        for c in self.coefficients.iter() {
            value = value + *c * x_i;
            x_i = x_i * point;
        }
        value
    }

    pub fn evaluate_domain(&self, domain: Vec<S>) -> Vec<S> {
        domain.into_iter().map(|p| self.evaluate(p)).collect()
    }

    /// Use lagrange interpolation to create a polynomial from lists of X and Y coordinates
    pub fn interpolate_domain(domain: Vec<S>, values: Vec<S>) -> Result<Self, InterpolationError> {
        if domain.len() != values.len() {
            return Err(InterpolationError::MismatchedCoordinateLists);
        }
        if domain.len() == 0 {
            return Err(InterpolationError::NoPointsProvided);
        }

        let x = Self::from(vec![S::zero(), S::one()]);
        let mut acc = Self::zero();
        for (i, v_i) in values.into_iter().enumerate() {
            let mut product = Polynomial::from(vec![v_i]);
            for (j, d_j) in domain.iter().enumerate() {
                dbg!(i, j);
                if j == i {
                    continue;
                }
                dbg!(&product);
                product = product
                    * ((x.clone() - Polynomial::from(vec![*d_j]))
                        * Polynomial::from(vec![(domain[i] - *d_j).inverse()]));
            }
            acc += product;
            dbg!(&acc);
        }
        Ok(acc)
    }

    pub fn zeroifier_domain(domain: Vec<S>) -> Self {
        let x = Self::from(vec![S::zero(), S::one()]);
        let mut acc = Self::from(vec![S::one()]);
        for d in domain {
            acc = acc * (x.clone() - Polynomial::from(vec![d]));
        }
        return acc;
    }

    pub fn scale(&self, factor: S) -> Self {
        Self::from(
            self.coefficients
                .iter()
                .enumerate()
                .map(|(i, coef)| factor.clone().pow(i) * coef.clone())
                .collect(),
        )
    }

    pub fn are_colinear(points: Vec<(S, S)>) -> bool {
        let mut domain = Vec::with_capacity(points.len());
        let mut range = Vec::with_capacity(points.len());
        for (x, y) in points {
            domain.push(x);
            range.push(y);
        }
        Self::interpolate_domain(domain, range)
            .expect("Provided the same number of x and y coordinates")
            .degree()
            == 1
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
            coefficients: self.coefficients.iter().map(|c| -*c).collect(),
        }
    }
}

impl<S> std::ops::Add for Polynomial<S>
where
    S: CopyScalar,
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

impl<S> std::ops::AddAssign for Polynomial<S>
where
    S: CopyScalar,
{
    fn add_assign(&mut self, mut rhs: Self) {
        if self.is_zero() {
            *self = rhs;
            return;
        } else if rhs.is_zero() {
            return;
        }

        if rhs.degree() > self.degree() {
            std::mem::swap(&mut self.coefficients, &mut rhs.coefficients)
        }
        for (c1, c2) in self.coefficients.iter_mut().zip(rhs.coefficients) {
            // TODO: adapt for AddAssign if it is included in Scalar
            *c1 = *c1 + c2
        }
    }
}

impl<S: CopyScalar> PartialOrd for Polynomial<S> {
    fn partial_cmp(&self, rhs: &Self) -> Option<Ordering> {
        if self.degree() > rhs.degree() {
            return Some(Ordering::Greater);
        }
        if self.degree() < rhs.degree() {
            return Some(Ordering::Less);
        }
        if self.is_zero() {
            return Some(Ordering::Equal);
        }
        for (self_coef, other_coef) in self
            .coefficients
            .iter()
            .rev()
            .zip(rhs.coefficients.iter().rev())
        {
            match self_coef.partial_cmp(other_coef) {
                Some(order) => match order {
                    Ordering::Less => return Some(Ordering::Less),
                    Ordering::Equal => continue,
                    Ordering::Greater => return Some(Ordering::Greater),
                },
                None => unreachable!(),
            }
        }
        Some(Ordering::Equal)
    }
}

impl<S> std::ops::Sub for Polynomial<S>
where
    S: CopyScalar,
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
    S: CopyScalar,
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
                // TODO: Adapt for AddAssign if it's added back
                coefs[i + j] = coefs[i + j] + *c1 * *c2
            }
        }

        return Polynomial::from(coefs);
    }
}

impl<S> std::ops::Mul<&Self> for Polynomial<S>
where
    S: CopyScalar,
{
    type Output = Self;

    fn mul(self, rhs: &Self) -> Self::Output {
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
                // TODO: Adapt for AddAssign if it's added back
                coefs[i + j] = coefs[i + j] + *c1 * *c2
            }
        }

        return Polynomial::from(coefs);
    }
}

// TODO: for better performance, we can use FFTs instead of
// polynomial long divison. For now, this will do.
impl<S> std::ops::Div for Polynomial<S>
where
    S: CopyScalar,
{
    type Output = Self;
    fn div(self, rhs: Self) -> Self::Output {
        let (quotient, _rem) = self.long_divide(rhs).expect("Cannot divide by zero");
        quotient
    }
}

impl<S> std::ops::Rem for Polynomial<S>
where
    S: CopyScalar,
{
    type Output = Self;

    fn rem(self, rhs: Self) -> Self::Output {
        let (_quotient, rem) = self.long_divide(rhs).expect("Cannot divide by zero");
        rem
    }
}

#[cfg(test)]
mod tests {

    use crate::{
        scalar::{CopyScalar, SimpleNum},
        GoldilocksElement,
    };

    impl CopyScalar for u32 {}
    use super::Polynomial;

    #[test]
    fn test_add() {
        let a = Polynomial::from_ref(&[1, 2, 3, 4, 5]);
        let b = Polynomial::from_ref(&[1, 2, 3]);

        assert_eq!(
            (a + b).coefficients,
            Polynomial::from_ref(&[2, 4, 6, 4, 5]).coefficients
        )
    }

    #[test]
    fn test_mul() {
        let a = Polynomial::from_ref(&[1, 2, 3, 4, 5]);
        let b = Polynomial::from_ref(&[1, 2, 3]);

        // [1, 2, 3, 4, 5]
        // [0, 2, 4, 6, 8, 10]
        // [0, 0, 3, 6, 9, 12, 15]
        assert_eq!(
            (a * b).coefficients,
            Polynomial::from_ref(&[1, 4, 10, 16, 22, 22, 15]).coefficients
        )
    }

    #[test]
    fn test_interpolate() {
        let domain: Vec<u32> = vec![0, 1, 2, 5];
        let values = vec![2u32, 3, 12, 147];

        let domain: Vec<GoldilocksElement> = domain.into_iter().map(|v| v.into()).collect();
        let values: Vec<GoldilocksElement> = values.into_iter().map(|v| v.into()).collect();

        let p =
            Polynomial::interpolate_domain(domain, values).expect("Polynomial must interpolate");
        p.evaluate(0.into());
        dbg!(&p);
        assert_eq!(p.evaluate(35.into()), 44067.into())
    }

    #[test]
    // fn test_interpolate_again() {
    //     let domain: Vec<u32> = vec![0, 1, 2];
    //     let values = vec![1u32, 2, 5];

    //     let domain: Vec<FieldElement<BigInt, DefaultField>> = domain
    //         .into_iter()
    //         .map(|v| FieldElement::new(v.into()))
    //         .collect();
    //     let values: Vec<FieldElement<BigInt, DefaultField>> = values
    //         .into_iter()
    //         .map(|v| FieldElement::new(v.into()))
    //         .collect();

    //     let p =
    //         Polynomial::interpolate_domain(domain, values).expect("Polynomial must interpolate");
    //     dbg!(&p);
    //     assert_eq!(
    //         p.evaluate(<FieldElement<BigInt, DefaultField>>::new(35u32.into())),
    //         <FieldElement<BigInt, DefaultField>>::new(1226.into())
    //     )
    // }
    #[test]
    fn test_scale() {
        let a = Polynomial::from(vec![1, 2, 3]);
        let c = a.scale(2).evaluate(10);
        assert_eq!(a.evaluate(20), c)
    }

    #[test]
    fn test_degree() {
        let a = Polynomial::from_ref(&[1, 2, 3]);
        let b: Polynomial<u32> = Polynomial::zero();
        assert_eq!(a.degree(), 2);
        assert_eq!(b.degree(), -1);
    }
}
