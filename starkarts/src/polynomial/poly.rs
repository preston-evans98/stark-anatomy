use std::cmp::Ordering;

use crate::{field::DivByZeroError, field_elem::FieldElement};

#[derive(Debug)]
pub enum InterpolationError {
    MismatchedCoordinateLists,
    NoPointsProvided,
}

#[derive(Debug, Clone, PartialEq)]
/// A polynomial, represented as a list of coefficients, with the lowest degree coefficient first.
/// For example, 5x^2 + 1 is represented [1, 0, 5]
pub struct Polynomial<F> {
    pub coefficients: Vec<F>,
}

impl<F> Polynomial<F>
where
    F: FieldElement,
{
    pub fn from<T: Into<F>>(coefs: Vec<T>) -> Self {
        let mut p = Self {
            coefficients: coefs.into_iter().map(|c| c.into()).collect(),
        };
        p._trim_degree();
        p
    }

    pub fn from_ref<T: Copy + Into<F>>(coefs: &[T]) -> Self {
        let mut p = Self {
            coefficients: coefs.iter().map(|c| (*c).into()).collect(),
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
    pub fn leading_coefficient(&self) -> Option<&F> {
        self.coefficients.last()
    }

    fn is_zero(&self) -> bool {
        self.coefficients.len() == 0
    }

    fn zero() -> Self {
        Self {
            coefficients: Vec::new(),
        }
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
        let mut quotient_coefs = F::zeros(result_degree_upper_bound);
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
                std::iter::repeat(F::zero())
                    .take(shift)
                    .chain(std::iter::once(coefficient.clone()))
                    .collect(),
            ) * divisor.clone();
            quotient_coefs[shift] = coefficient;
            remainder = remainder - subtrahend;
        }
        let quotient = Polynomial::from(quotient_coefs);
        return Ok((quotient, remainder));
    }

    /// A naive implementation of evaluation at a point
    pub fn evaluate<T: Into<F>>(&self, point: T) -> F {
        let point = point.into();
        let mut x_i = F::one();
        let mut value = F::zero();
        for c in self.coefficients.iter() {
            value = value + c.clone() * x_i.clone();
            x_i = x_i * point;
        }
        value
    }

    pub fn evaluate_domain<T: Into<F>>(&self, domain: Vec<T>) -> Vec<F> {
        domain.into_iter().map(|p| self.evaluate(p)).collect()
    }

    /// Use lagrange interpolation to create a polynomial from lists of X and Y coordinates
    pub fn interpolate_domain<T: Into<F>>(
        domain: Vec<T>,
        values: Vec<T>,
    ) -> Result<Self, InterpolationError> {
        if domain.len() != values.len() {
            return Err(InterpolationError::MismatchedCoordinateLists);
        }
        if domain.len() == 0 {
            return Err(InterpolationError::NoPointsProvided);
        }

        let values: Vec<F> = values.into_iter().map(|v| v.into()).collect();
        let domain: Vec<F> = domain.into_iter().map(|v| v.into()).collect();

        let x = Self::from(vec![F::zero(), F::one()]);
        let mut acc = Self::zero();
        for (i, v_i) in values.into_iter().enumerate() {
            let mut product = Polynomial::from(vec![v_i]);
            for (j, d_j) in domain.iter().enumerate() {
                if j == i {
                    continue;
                }
                product = product
                    * (x.clone() - Polynomial::from(vec![*d_j]))
                    * Polynomial::from(vec![(domain[i] - *d_j).inverse()]);
            }
            acc += product;
            dbg!(&acc);
        }
        Ok(acc)
    }

    pub fn zeroifier_domain(domain: Vec<F>) -> Self {
        let x = Self::from(vec![F::zero(), F::one()]);
        let mut acc = Self::from(vec![F::one()]);
        for d in domain {
            acc = acc * (x.clone() - Polynomial::from(vec![d]));
        }
        return acc;
    }

    pub fn scale<T: Into<F>>(&self, factor: T) -> Self {
        let factor = factor.into();
        Self::from(
            self.coefficients
                .iter()
                .enumerate()
                .map(|(i, coef)| factor.pow(i) * (*coef))
                .collect(),
        )
    }

    pub fn are_colinear(points: Vec<(F, F)>) -> bool {
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

    // Trim leading zero coefficients from the polynomial, so that future calls
    // to degree() return the correct result
    fn _trim_degree(&mut self) {
        while self.coefficients.ends_with(&[F::zero()]) {
            self.coefficients.pop();
        }
    }
}

impl<F> std::ops::Neg for Polynomial<F>
where
    F: FieldElement,
{
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self {
            coefficients: self.coefficients.iter().map(|c| -*c).collect(),
        }
    }
}

impl<F> std::ops::Add for Polynomial<F>
where
    F: FieldElement,
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
                .zip(short_iter.chain(std::iter::repeat(F::zero())))
                .map(|(l, r)| l + r)
                .collect(),
        );
    }
}

impl<F> std::ops::AddAssign for Polynomial<F>
where
    F: FieldElement,
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
            // TODO: adapt for AddAssign if it is included in FieldElement
            *c1 = c1.clone() + c2
        }
    }
}

impl<F: FieldElement> PartialOrd for Polynomial<F> {
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

impl<F> std::ops::Sub for Polynomial<F>
where
    F: FieldElement,
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
                .zip(short_iter.chain(std::iter::repeat(F::zero())))
                .map(|(l, r)| l - r)
                .collect(),
        );
    }
}

impl<F> std::ops::Mul for Polynomial<F>
where
    F: FieldElement,
{
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        if self.is_zero() || rhs.is_zero() {
            return Self::zero();
        }
        let mut coefs = F::zeros((self.degree() + rhs.degree() + 1) as usize);

        for (i, c1) in self.coefficients.iter().enumerate() {
            // If c1 is zero then c1 * c2 = 0 * c2 = 0, so we can skip this iteration
            if c1.is_zero() {
                continue;
            }
            for (j, c2) in rhs.coefficients.iter().enumerate() {
                // TODO: Adapt for AddAssign if it's added back
                coefs[i + j] = coefs[i + j].clone() + c1.clone() * c2.clone()
            }
        }

        return Polynomial::from(coefs);
    }
}

// TODO: for better performance, we can use FFTs instead of
// polynomial long divison. For now, this will do.
impl<F> std::ops::Div for Polynomial<F>
where
    F: FieldElement,
{
    type Output = Self;
    fn div(self, rhs: Self) -> Self::Output {
        let (quotient, _rem) = self.long_divide(rhs).expect("Cannot divide by zero");
        quotient
    }
}

impl<F> std::ops::Rem for Polynomial<F>
where
    F: FieldElement,
{
    type Output = Self;

    fn rem(self, rhs: Self) -> Self::Output {
        let (_quotient, rem) = self.long_divide(rhs).expect("Cannot divide by zero");
        rem
    }
}

#[cfg(test)]
mod tests {

    use crate::{field_elem::FieldElement, GoldilocksElement};

    use super::Polynomial;

    impl<F: FieldElement, T: Into<F> + Copy> PartialEq<Vec<T>> for Polynomial<F> {
        fn eq(&self, other: &Vec<T>) -> bool {
            if self.coefficients.len() != other.len() {
                return false;
            }
            for (c1, c2) in self.coefficients.iter().zip(other.iter()) {
                if *c1 != (*c2).into() {
                    return false;
                }
            }
            true
        }
    }

    #[test]
    fn test_add() {
        let a: Polynomial<GoldilocksElement> = Polynomial::from_ref(&[1, 2, 3, 4, 5]);
        let b = Polynomial::from_ref(&[1, 2, 3]);

        assert_eq!((a + b), vec![2, 4, 6, 4, 5])
    }

    #[test]
    fn test_mul() {
        let a: Polynomial<GoldilocksElement> = Polynomial::from_ref(&[1, 2, 3, 4, 5]);
        let b = Polynomial::from_ref(&[1, 2, 3]);

        // [1, 2, 3, 4, 5]
        // [0, 2, 4, 6, 8, 10]
        // [0, 0, 3, 6, 9, 12, 15]
        assert_eq!(a * b, vec![1, 4, 10, 16, 22, 22, 15])
    }

    #[test]
    fn test_interpolate() {
        let domain: Vec<u32> = vec![0, 1, 2, 5];
        let values = vec![2u32, 3, 12, 147];

        let domain: Vec<GoldilocksElement> = domain.iter().map(|v| (*v).into()).collect();
        let values: Vec<GoldilocksElement> = values.iter().map(|v| (*v).into()).collect();

        let p: Polynomial<GoldilocksElement> =
            Polynomial::interpolate_domain(domain, values).expect("Polynomial must interpolate");
        dbg!(&p);
        assert_eq!(p.evaluate(35), 44067)
    }
    #[test]
    fn test_interpolate_again() {
        let domain: Vec<u32> = vec![0, 1, 2];
        let values = vec![1u32, 2, 5];

        let domain: Vec<GoldilocksElement> = domain.iter().map(|v| (*v).into()).collect();
        let values: Vec<GoldilocksElement> = values.iter().map(|v| (*v).into()).collect();

        let p: Polynomial<GoldilocksElement> =
            Polynomial::interpolate_domain(domain, values).expect("Polynomial must interpolate");
        dbg!(&p);
        assert_eq!(p.evaluate(35), 1226)
    }
    // #[test]
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
        let a: Polynomial<GoldilocksElement> = Polynomial::from(vec![1, 2, 3]);
        let c = a.scale(2).evaluate(10);
        assert_eq!(a.scale(2), vec![1, 4, 12]);
        assert_eq!(a.evaluate(20), c)
    }

    #[test]
    fn test_degree() {
        let a: Polynomial<GoldilocksElement> = Polynomial::from_ref(&[1u32, 2, 3]);
        let b: Polynomial<GoldilocksElement> = Polynomial::zero();
        assert_eq!(a.degree(), 2);
        assert_eq!(b.degree(), -1);
    }
}
