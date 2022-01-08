use std::ops::{Add, AddAssign, Div, Mul, Neg, Rem, Sub};

use crate::{scalar::Scalar, PrimeField};

#[derive(Debug, Clone)]
pub struct FieldElement<S, F> {
    value: S,
    field: std::marker::PhantomData<F>,
}

impl<S, F> Scalar for FieldElement<S, F>
where
    S: Scalar,
    F: PrimeField<Elem = S>,
{
    fn is_zero(&self) -> bool {
        self.value.is_zero()
    }

    fn zero() -> Self {
        FieldElement::new(S::zero())
    }

    fn one() -> Self {
        FieldElement::new(S::one())
    }
}

impl<S, F> Rem for FieldElement<S, F>
where
    S: Scalar,
    F: PrimeField<Elem = S>,
{
    type Output = Self;

    fn rem(self, rhs: Self) -> Self::Output {
        FieldElement {
            value: self.value % rhs.value,
            field: std::marker::PhantomData,
        }
    }
}

impl<S, F> Sub for FieldElement<S, F>
where
    S: Scalar,
    F: PrimeField<Elem = S>,
{
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        FieldElement::new(F::subtract(self.value, rhs.value))
    }
}

impl<S, F> Neg for FieldElement<S, F>
where
    S: Scalar,
    F: PrimeField<Elem = S>,
{
    type Output = Self;
    fn neg(self) -> Self::Output {
        FieldElement::new(F::subtract(F::zero(), self.value))
    }
}

impl<S, F> Mul for FieldElement<S, F>
where
    S: Scalar,
    F: PrimeField<Elem = S>,
{
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        FieldElement::new(F::multiply(self.value, rhs.value))
    }
}

impl<S, F> Div for FieldElement<S, F>
where
    S: Scalar,
    F: PrimeField<Elem = S>,
{
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        FieldElement::new(F::divide(self.value, rhs.value).expect("Cannot divide by zero"))
    }
}

impl<S, F> PartialEq for FieldElement<S, F>
where
    S: Scalar,
    F: PrimeField<Elem = S>,
{
    fn eq(&self, other: &Self) -> bool {
        self.value == other.value
    }
}

impl<S, F> PartialOrd for FieldElement<S, F>
where
    S: Scalar,
    F: PrimeField<Elem = S>,
{
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.value.partial_cmp(&other.value)
    }
}

impl<S, F> FieldElement<S, F>
where
    S: Scalar,
    F: PrimeField<Elem = S>,
{
    pub fn new(value: S) -> Self {
        Self {
            value,
            field: std::marker::PhantomData,
        }
    }

    pub fn pow(&self, exponent: usize) -> Self {
        Self::pow_iter(self.value.clone(), F::one(), exponent)
    }

    fn pow_iter(value: S, acc: S, exponent: usize) -> Self {
        match exponent {
            0 => FieldElement {
                value: F::one(),
                field: std::marker::PhantomData,
            },
            1 => FieldElement {
                value: F::multiply(value, acc),
                field: std::marker::PhantomData,
            },
            x if x & 1 == 0 => {
                Self::pow_iter(F::multiply(value.clone(), value), acc, exponent >> 1)
            }
            x if x & 1 == 1 => Self::pow_iter(value.clone(), F::multiply(acc, value), exponent - 1),
            _ => unreachable!(),
        }
    }
}

impl<S, F> AddAssign<S> for FieldElement<S, F>
where
    S: Scalar,
    F: PrimeField<Elem = S>,
{
    fn add_assign(&mut self, rhs: S) {
        self.value = F::add(self.value.clone(), rhs)
    }
}

impl<S, F> AddAssign<FieldElement<S, F>> for FieldElement<S, F>
where
    S: Scalar,
    F: PrimeField<Elem = S>,
{
    fn add_assign(&mut self, rhs: FieldElement<S, F>) {
        self.value = F::add(self.value.clone(), rhs.value.clone())
    }
}

impl<S, F> Add<S> for FieldElement<S, F>
where
    S: Scalar,
    F: PrimeField<Elem = S>,
{
    type Output = FieldElement<S, F>;
    fn add(self, rhs: S) -> Self::Output {
        let value = F::add(self.value.clone(), rhs);
        FieldElement {
            value,
            field: self.field,
        }
    }
}

impl<S, F> Add<Self> for FieldElement<S, F>
where
    S: Scalar,
    F: PrimeField<Elem = S>,
{
    type Output = FieldElement<S, F>;
    fn add(self, rhs: Self) -> Self::Output {
        let value = F::add(self.value.clone(), rhs.value.clone());
        FieldElement {
            value,
            field: self.field,
        }
    }
}

#[cfg(test)]
mod tests {
    use num::{BigUint, Num};

    use crate::DefaultField;

    use super::FieldElement;

    #[test]
    fn pow_test() {
        let elem: FieldElement<BigUint, DefaultField> = FieldElement::new(12345u16.into());
        let result = elem.pow(2345);
        assert_eq!(
            result.value,
            BigUint::from_str_radix("23520667101221308275390517182423299422", 10).unwrap()
        )
    }
}
