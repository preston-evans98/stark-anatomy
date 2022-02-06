use std::{collections::HashMap, hash::Hash, mem::MaybeUninit};

use crate::{field_elem::FieldElement, Polynomial};

#[derive(Debug, Copy, Clone, Hash, PartialEq, Eq)]
pub struct Exponents<const VARIABLES: usize>([usize; VARIABLES]);

impl<const VARIABLES: usize> Exponents<VARIABLES> {
    const fn new() -> Exponents<VARIABLES> {
        Self([0usize; VARIABLES])
    }

    fn with_val_at_idx(val: usize, idx: usize) -> Self {
        let mut res = [0usize; VARIABLES];
        res[idx] = val;
        Self(res)
    }
}

impl<const VARIABLES: usize> std::ops::Mul for Exponents<VARIABLES> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let mut result = self;
        for (res, r) in result.0.iter_mut().zip(rhs.0.iter()) {
            *res += *r
        }
        result
    }
}

#[derive(Debug, Clone)]
/// A polynomial, represented as a set of tuples , with the lowest degree coefficient first.
/// For example, 5x^2 + 1 is represented [1, 0, 5]
pub struct MultiVariatePoly<F, const VARIABLES: usize> {
    dict: HashMap<Exponents<VARIABLES>, F>,
}

impl<F: FieldElement, const VARIABLES: usize> MultiVariatePoly<F, VARIABLES> {
    /// Returns the zero value of this polynomial's exponent representation
    fn zero_exponent() -> Exponents<VARIABLES> {
        Exponents::<VARIABLES>::new()
    }

    /// Returns a constant as a polynomial
    pub fn constant(c: F) -> Self {
        Self {
            dict: HashMap::from([(Self::zero_exponent(), c)]),
        }
    }

    /// Returns a polynomial with a single term
    pub fn monomial(exp: Exponents<VARIABLES>, coef: F) -> Self {
        Self {
            dict: HashMap::from([(exp, coef)]),
        }
    }

    /// Returns the zero polynomial (additive identity)
    pub fn zero() -> Self {
        Self {
            dict: HashMap::new(),
        }
    }
    /// Create a polynomial whose underlying storage has capacity for `capacity` terms
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            dict: HashMap::with_capacity(capacity),
        }
    }

    /// Raises this polynomial to a power
    pub fn pow(&self, mut exponent: usize) -> Self {
        if self.is_zero() {
            return Self::zero();
        }
        let mut acc = MultiVariatePoly::constant(F::one());
        let mut base = self.clone();
        while exponent != 0 {
            if exponent & 1 == 0 {
                base = &base * &base;
                exponent >>= 1
            } else {
                acc = &base * &acc;
                exponent -= 1
            }
        }
        acc
    }

    /// Interprets a univariate polynomial as multivariate. The caller must specify
    /// the index of the univariate's variable into the multivariate's exponent array.
    pub fn lift(polynomial: Polynomial<F>, index_into_exponents: usize) -> Self {
        let mut result = Self::with_capacity((polynomial.degree() + 1) as usize);
        for (power, coef) in polynomial.coefficients.iter().enumerate() {
            let exponents = Exponents::with_val_at_idx(power, index_into_exponents);
            result.insert(exponents, *coef);
        }
        result
    }

    /// Evaluates this polynomial at the provided points
    pub fn evaluate(&self, at: [F; VARIABLES]) -> F {
        let mut acc = F::zero();
        for (exp, coef) in self.dict.iter() {
            acc = acc
                + at.into_iter()
                    .zip(exp.0.into_iter())
                    .fold(*coef, |accumulator, (point, exponent)| {
                        accumulator * point.pow(exponent)
                    });
        }
        acc
    }

    pub fn evaluate_symbolic(&self, at: [Polynomial<F>; VARIABLES]) -> Polynomial<F> {
        let mut acc = Polynomial::zero();
        for (exponents, coef) in self.dict.iter() {
            let mut prod = Polynomial::<F>::from_ref(&[*coef]);
            for (exp, point) in exponents.0.iter().zip(at.iter()) {
                prod = prod * point.pow(*exp);
            }
            acc = acc + prod
        }
        acc
    }
}

impl<F: FieldElement, const VARIABLES: usize> MultiVariatePoly<F, VARIABLES> {
    /// Insert a term into the polynomial, overwriting whatever coefficient it might have had before.
    fn insert(&mut self, exponents: Exponents<VARIABLES>, value: F) {
        if value.is_zero() {
            self.dict.remove(&exponents);
        }
        self.dict.insert(exponents, value);
    }

    /// Add a single term to the polynomial. For example (3xy + 4x).add_term(2x) = (3xy + 6x)
    pub fn add_term(&mut self, exponents: Exponents<VARIABLES>, coef: F) {
        if coef.is_zero() {
            return;
        }
        match self.dict.insert(exponents, coef) {
            Some(existing_coef) => {
                let new_coef = existing_coef + coef;
                if new_coef.is_zero() {
                    self.dict.remove(&exponents);
                    return;
                }
                self.dict.insert(exponents, new_coef);
            }
            None => {}
        }
    }

    pub fn is_zero(&self) -> bool {
        self.len() == 0
    }

    /// The number of non-zero terms in the polynomial
    pub fn len(&self) -> usize {
        self.dict.len()
    }

    /// Return an array of monomials, each consisting of a single variable from the polynomial
    /// For example, the multivariate polynomial class [c x^a y^b z^c] returns [x, y, z]
    pub fn variables() -> [Self; VARIABLES] {
        let mut arr: [MaybeUninit<Self>; VARIABLES] =
            unsafe { MaybeUninit::uninit().assume_init() };
        let mut idx = 0;
        for elem in &mut arr {
            let exp = Exponents::<VARIABLES>::with_val_at_idx(1, idx);
            elem.write(Self::monomial(exp, F::one()));
            idx += 1;
        }

        // Hack for transmuting const generic arrays taken from: https://github.com/rust-lang/rust/issues/61956
        unsafe {
            let ptr = &mut arr as *mut _ as *mut [Self; VARIABLES];
            let res = ptr.read();
            core::mem::forget(arr);
            res
        }
    }
}

impl<F, const VARIABLES: usize> std::ops::Add for &MultiVariatePoly<F, VARIABLES>
where
    F: FieldElement,
{
    type Output = MultiVariatePoly<F, VARIABLES>;

    fn add(self, rhs: Self) -> Self::Output {
        let mut result = MultiVariatePoly::with_capacity(self.len() + rhs.len());
        for (exp, coef) in self.dict.iter() {
            result.insert(*exp, *coef);
        }
        for (exp, coef) in rhs.dict.iter() {
            result.add_term(*exp, *coef);
        }
        result
    }
}

impl<F, const VARIABLES: usize> std::ops::Mul for &MultiVariatePoly<F, VARIABLES>
where
    F: FieldElement,
{
    type Output = MultiVariatePoly<F, VARIABLES>;

    fn mul(self, rhs: Self) -> Self::Output {
        let mut result = MultiVariatePoly::with_capacity(self.len() * rhs.len());
        for (exp1, coef1) in self.dict.iter() {
            for (exp2, coef2) in rhs.dict.iter() {
                let coef = *coef1 * *coef2;
                let exp = *exp1 * *exp2;
                result.add_term(exp, coef)
            }
        }
        result
    }
}

impl<F: FieldElement, const VARIABLES: usize> std::ops::Neg for &MultiVariatePoly<F, VARIABLES> {
    type Output = MultiVariatePoly<F, VARIABLES>;

    fn neg(self) -> Self::Output {
        let mut other_dict = HashMap::with_capacity(self.len());
        for (k, v) in self.dict.iter() {
            other_dict.insert(*k, -*v);
        }
        MultiVariatePoly { dict: other_dict }
    }
}

impl<F, const VARIABLES: usize> std::ops::Sub for &MultiVariatePoly<F, VARIABLES>
where
    F: FieldElement,
{
    type Output = MultiVariatePoly<F, VARIABLES>;

    fn sub(self, rhs: Self) -> Self::Output {
        let mut result = MultiVariatePoly::with_capacity(self.len() + rhs.len());
        for (exp, coef) in self.dict.iter() {
            result.insert(*exp, -*coef);
        }
        for (exp, coef) in rhs.dict.iter() {
            result.add_term(*exp, -*coef)
        }
        result
    }
}
