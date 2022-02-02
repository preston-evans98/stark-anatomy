use std::collections::HashMap;

use crate::{field_elem::FieldElement, P};
type exponent = usize;

#[derive(Debug, Copy, Clone, Hash, PartialEq, Eq)]
pub struct Exponents<const VARIABLES: usize>([exponent; VARIABLES]);

impl<const VARIABLES: usize> Exponents<VARIABLES> {
    const fn new() -> Exponents<VARIABLES> {
        Self([0usize; VARIABLES])
    }
    fn pad<const NEW_LEN: usize>(&self) -> Exponents<NEW_LEN> {
        let mut result = [0; NEW_LEN];
        result[..VARIABLES].copy_from_slice(&self.0);
        Exponents(result)
    }
    const fn get_size() -> usize {
        VARIABLES
    }
}

impl<const VARIABLES: usize, const OTHER: usize> std::ops::Mul<Exponents<OTHER>>
    for Exponents<VARIABLES>
where
    [(); max(VARIABLES, OTHER)]:,
{
    type Output = Exponents<{ max(VARIABLES, OTHER) }>;

    fn mul(self, rhs: Exponents<OTHER>) -> Self::Output {
        let mut result = <Self as std::ops::Mul<Exponents<OTHER>>>::Output::new();
        for (res, l) in result.0.iter_mut().zip(self.0.iter()) {
            *res += *l
        }
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

impl<F, const VARIABLES: usize> MultiVariatePoly<F, VARIABLES> {
    fn insert_padded<const OTHER: usize>(&mut self, exponents: Exponents<OTHER>, value: F) {
        let exponents = exponents.pad::<VARIABLES>();
        self.dict.insert(exponents, value);
    }

    unsafe fn unsafe_insert<const OTHER: usize>(
        &mut self,
        mut exponents: Exponents<OTHER>,
        value: F,
    ) -> Option<F> {
        // Hack for tranmuting cons generic arrays taken from: https://github.com/rust-lang/rust/issues/61956
        let ptr = &mut exponents.0 as *mut _ as *mut [exponent; VARIABLES];
        let res = ptr.read();
        self.dict.insert(Exponents(res), value)
    }

    // fn insert_le(&mut self, exponents: Exponents<VARIABLES>, value: F) {
    //     if VARIABLES < OTHER {
    //         let padded = exponents.pad::<OTHER>();
    //     }
    //     self.dict.insert(exponents, value);
    // }
}

impl<F, const VARIABLES: usize> MultiVariatePoly<F, VARIABLES> {
    pub fn zero() -> Self {
        Self {
            dict: HashMap::new(),
        }
    }
}
pub const fn max(lhs: usize, rhs: usize) -> usize {
    if lhs > rhs {
        lhs
    } else {
        rhs
    }
}

impl<F, const VARIABLES: usize, const OTHER: usize> std::ops::Add<&MultiVariatePoly<F, OTHER>>
    for &MultiVariatePoly<F, VARIABLES>
where
    [(); max(VARIABLES, OTHER)]:,
    F: FieldElement,
{
    type Output = MultiVariatePoly<F, { max(VARIABLES, OTHER) }>;

    fn add(self, rhs: &MultiVariatePoly<F, OTHER>) -> Self::Output {
        let mut result = Self::Output::zero();
        if VARIABLES > OTHER {
            for (exponents, coefficient) in rhs.dict.iter() {
                result.insert_padded(*exponents, *coefficient);
            }
            for (exponents, coefficient) in self.dict.iter() {
                unsafe {
                    match result.unsafe_insert(*exponents, *coefficient) {
                        Some(existing) => {
                            result.unsafe_insert(*exponents, *coefficient + existing);
                        }
                        None => {}
                    }
                }
            }
            return result;
        } else if VARIABLES == OTHER {
            for (exponents, coefficient) in rhs.dict.iter() {
                unsafe {
                    result.unsafe_insert(*exponents, *coefficient);
                }
            }
            for (exponents, coefficient) in self.dict.iter() {
                unsafe {
                    match result.unsafe_insert(*exponents, *coefficient) {
                        Some(existing) => {
                            result.unsafe_insert(*exponents, *coefficient + existing);
                        }
                        None => {}
                    }
                }
            }
        } else {
            for (exponents, coefficient) in self.dict.iter() {
                result.insert_padded(*exponents, *coefficient)
            }
            for (exponents, coefficient) in rhs.dict.iter() {
                unsafe {
                    match result.unsafe_insert(*exponents, *coefficient) {
                        Some(existing) => {
                            result.unsafe_insert(*exponents, *coefficient + existing);
                        }
                        None => {}
                    }
                }
            }
        }

        result
    }
}

// impl<F> std::iter::Iterator for MultiVariatePoly<F> {
//     type Item = (Variable, Vec<F>);

//     fn next(&mut self) -> Option<Self::Item> {

//     }
// }

impl<F, const VARIABLES: usize, const OTHER: usize> std::ops::Mul<&MultiVariatePoly<F, OTHER>>
    for &MultiVariatePoly<F, VARIABLES>
where
    [(); max(VARIABLES, OTHER)]:,
    F: FieldElement,
{
    type Output = MultiVariatePoly<F, { max(VARIABLES, OTHER) }>;

    fn mul(self, rhs: &MultiVariatePoly<F, OTHER>) -> Self::Output {
        let mut result = MultiVariatePoly::zero();
        for (exp1, coef1) in self.dict.iter() {
            for (exp2, coef2) in rhs.dict.iter() {
                let coef = *coef1 * *coef2;
                let exp = *exp1 * *exp2;
                match result.dict.insert(exp, coef) {
                    Some(existing_coef) => {
                        result.dict.insert(exp, existing_coef + coef);
                    }
                    None => {}
                }
            }
        }
        result
    }
}
