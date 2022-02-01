use std::collections::HashMap;

use crate::field_elem::FieldElement;
type exponent = usize;

#[derive(Debug, Copy, Clone, Hash, PartialEq, Eq)]
pub struct Exponents<const NUM_VARIABLES: usize>([exponent; NUM_VARIABLES]);
// impl<const NUM_VARIABLES: usize> Exponents<NUM_VARIABLES> {
//     fn pad_self<const TO_LEN: usize>(self) -> Exponents<TO_LEN> {
//         let mut result = [0; TO_LEN];
//         result[..NUM_VARIABLES].copy_from_slice(&self.0);
//         Exponents(result)
//     }
// }
pub trait Paddable<const OTHER: usize> {
    fn pad(smaller: Exponents<OTHER>) -> Self;
}

impl<const VARIABLES: usize, const OTHER: usize> Paddable<OTHER> for Exponents<VARIABLES>
where
    Compare<VARIABLES, OTHER>: Greater,
    Assert<{ VARIABLES > OTHER }>: IsTrue,
{
    fn pad(smaller: Exponents<OTHER>) -> Exponents<VARIABLES> {
        let mut result = [0; VARIABLES];
        result[..OTHER].copy_from_slice(&smaller.0);
        Exponents(result)
    }
}

impl<const VARIABLES: usize, const OTHER: usize> Paddable<OTHER> for Exponents<VARIABLES>
where
    IsLessThan<OTHER, VARIABLES>: IsTrue2,
{
    fn pad(smaller: Exponents<OTHER>) -> Exponents<VARIABLES> {
        let mut result = [0; VARIABLES];
        result[..OTHER].copy_from_slice(&smaller.0);
        Exponents(result)
    }
}

// impl<const VARIABLES: usize, const SMALLER: usize> std::ops::Add<Exponents<SMALLER>>
//     for Exponents<VARIABLES>
// {
//     type Output = Exponents<VARIABLES>;

//     fn add(self, rhs: Exponents<SMALLER>) -> Self::Output {
//         let mut result: [exponent; VARIABLES] = [0; VARIABLES];
//         for (result_elem, (lhs, rhs)) in result.iter_mut().zip(self.0.iter().zip(rhs.0.iter())) {
//             *result_elem = *lhs + *rhs;
//         }
//         Exponents(result)
//     }
// }

// impl<const VARIABLES: usize, const SMALLER: usize> std::ops::AddAssign<Exponents<SMALLER>>
//     for Exponents<VARIABLES>
// {
//     type Output = Exponents<VARIABLES>;

//     fn add(self, rhs: Exponents<SMALLER>) -> Self::Output {
//         let mut result: [exponent; VARIABLES] = [0; VARIABLES];
//         for (result_elem, (lhs, rhs)) in result.iter_mut().zip(self.0.iter().zip(rhs.0.iter())) {
//             *result_elem = *lhs + *rhs;
//         }
//         Exponents(result)
//     }
// }

#[derive(Debug, Clone)]
/// A polynomial, represented as a set of tuples , with the lowest degree coefficient first.
/// For example, 5x^2 + 1 is represented [1, 0, 5]
pub struct MultiVariatePoly<F, const VARIABLES: usize> {
    dict: HashMap<Exponents<VARIABLES>, F>,
}

pub struct Compare<const LEFT: usize, const RIGHT: usize> {}
pub trait Greater {}
impl<const LEFT: usize, const RIGHT: usize> Greater for Compare<LEFT, RIGHT> where
    Assert<{ LEFT > RIGHT }>: IsTrue
{
}

pub trait Lesser {}
impl<const LEFT: usize, const RIGHT: usize> Lesser for Compare<LEFT, RIGHT> where
    Assert2<{ LEFT < RIGHT }>: IsTrue
{
}

pub trait Marker {}

pub trait Equal {}
impl<const LEFT: usize, const RIGHT: usize> Equal for Compare<LEFT, RIGHT> where
    Assert<{ LEFT == RIGHT }>: IsTrue
{
}

pub struct IsLessThan<const LEFT: usize, const RIGHT: usize> {}
impl<const LEFT: usize, const RIGHT: usize> IsTrue2 for IsLessThan<LEFT, RIGHT> where
    Assert<{ LEFT < RIGHT }>: IsTrue2
{
}
pub struct Assert<const COND: bool> {}
pub struct Assert2<const COND: bool> {}

pub trait IsTrue {}
pub trait IsTrue2 {}
pub trait IsFalse {}

impl IsTrue for Assert<true> {}
impl IsTrue2 for Assert2<true> {}
impl IsFalse for Assert<false> {}

// trait LessThan<o

impl<F, const VARIABLES: usize> MultiVariatePoly<F, VARIABLES>
// where
//     [(); max(VARIABLES, OTHER)]:,
{
    // Where other <= VARIABLES
    // fn insert_le<const OTHER: usize>(&mut self, exponents: Exponents<OTHER>, value: F) {
    //     let padded = exponents.pad::<VARIABLES>();
    //     self.dict.insert(padded, value);
    // }
    fn insert_padded<const SMALLER: usize>(
        &mut self,
        exponents: Exponents<SMALLER>,
        value: F,
    ) -> Option<F>
    where
        Compare<VARIABLES, SMALLER>: Greater,
        Assert<{ VARIABLES > SMALLER }>: IsTrue,
        // [(); { VARIABLES > OTHER }]:,
    {
        let exponents = Exponents::<VARIABLES>::pad(exponents);
        self.dict.insert(exponents, value)
    }

    fn insert_padded2<const SMALLER: usize>(
        &mut self,
        exponents: Exponents<SMALLER>,
        value: F,
    ) -> Option<F>
    where
        IsLessThan<SMALLER, VARIABLES>: IsTrue2, // [(); { VARIABLES > OTHER }]:,
        Assert2<{ SMALLER < VARIABLES }>: IsTrue2,
    {
        let exponents = Exponents::<VARIABLES>::pad(exponents);
        self.dict.insert(exponents, value)
    }

    // fn get_mut_padded<const SMALLER: usize>(
    //     &mut self,
    //     exponents: Exponents<SMALLER>,
    // ) -> Option<&mut F>
    // where
    //     Compare<VARIABLES, SMALLER>: Greater,
    //     // [(); { VARIABLES > OTHER }]:,
    // {
    //     let exponents = Exponents::<VARIABLES>::pad(exponents);
    //     self.dict.get_mut(&exponents)
    // }

    fn insert(&mut self, exponents: Exponents<VARIABLES>, value: F) -> Option<F> {
        self.dict.insert(exponents, value)
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

impl<F, const VARIABLES: usize, const SMALLER: usize> std::ops::Add<&MultiVariatePoly<F, SMALLER>>
    for &MultiVariatePoly<F, VARIABLES>
where
    Compare<VARIABLES, SMALLER>: Greater,
    Assert<{ VARIABLES > SMALLER }>: IsTrue,
    F: FieldElement,
{
    type Output = MultiVariatePoly<F, VARIABLES>;

    fn add(self, rhs: &MultiVariatePoly<F, SMALLER>) -> Self::Output {
        let mut result: MultiVariatePoly<F, VARIABLES> = MultiVariatePoly::zero();
        for (exponents, coefficient) in rhs.dict.iter() {
            result.insert_padded(*exponents, *coefficient);
        }
        for (exponents, coefficient) in self.dict.iter() {
            match result.insert(*exponents, *coefficient) {
                Some(item) => {
                    result.insert(*exponents, *coefficient + item);
                }
                None => {}
            };
        }
        result
    }
}

impl<F, const VARIABLES: usize, const BIGGER: usize> std::ops::Add<&MultiVariatePoly<F, BIGGER>>
    for &MultiVariatePoly<F, VARIABLES>
where
    IsLessThan<VARIABLES, BIGGER>: IsTrue2,
    Assert2<{ VARIABLES < BIGGER }>: IsTrue2,
    F: FieldElement,
{
    type Output = MultiVariatePoly<F, BIGGER>;

    fn add(self, rhs: &MultiVariatePoly<F, BIGGER>) -> Self::Output {
        let mut result: MultiVariatePoly<F, BIGGER> = MultiVariatePoly::zero();
        for (small_exponents, coefficient) in self.dict.iter() {
            match result.insert_padded2(*small_exponents, *coefficient) {
                Some(item) => {
                    result.insert_padded2(*small_exponents, *coefficient + item);
                }
                None => {}
            };
        }
        for (exponents, coefficient) in rhs.dict.iter() {
            result.insert(*exponents, *coefficient);
        }

        result
    }
}

// impl<F> std::iter::Iterator for MultiVariatePoly<F> {
//     type Item = (Variable, Vec<F>);

//     fn next(&mut self) -> Option<Self::Item> {

//     }
// }
