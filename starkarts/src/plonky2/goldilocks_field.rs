use crate::{field::PrimeField, uint::Uint};

#[derive(Debug, Clone, Copy)]
pub struct GoldilocksField;
pub const P: u64 = u64::MAX - (u32::MAX as u64) + 1;
const LOW_32_BITS: u64 = std::u32::MAX as u64;
const U32_MAX: u64 = std::u32::MAX as u64;
// const U32_MAX_DOUBLED: u64 = (std::u32::MAX as u64) * 2;

///  A 64 bit field allowing very fast arithmetic on modern CPUs.
///  https://github.com/mir-protocol/plonky2/blob/main/plonky2.pdf
impl GoldilocksField {
    /// A (somewhat) optimized function to compute X % p for any 128 bit integer x,
    /// Based on section 2 of https://github.com/mir-protocol/plonky2/blob/main/plonky2.pdf
    /// Note that this function computes a reduction to a u64, but does not put the result into
    /// canonical form. To get the element into canonical form, a single subtraction
    /// of P may be necessary.
    pub fn partial_reduce(low_order_bits: u64, high_order_bits: u64) -> u64 {
        let n0 = low_order_bits;
        let n1: u64 = high_order_bits & LOW_32_BITS;
        let n1: u64 = (n1 << 32) - n1;
        let n2 = high_order_bits >> 32;
        let (mut intermediate, overflow) = n0.overflowing_add(n1);
        if overflow {
            // If we overflowed, this subtraction will underflow yielding the correct result
            intermediate = Self::sub_p_from(intermediate)
        }
        let (mut result, underflow) = intermediate.overflowing_sub(n2);
        if underflow {
            result = result.wrapping_add(P)
        }
        result
    }

    pub fn reduce_multiplication(widening_mul_result: (u64, u64)) -> u64 {
        Self::partial_reduce(widening_mul_result.0, widening_mul_result.1)
    }
    /// Reduce a u128 where the highest 32 bits are zero. This function is slightly more optimized than
    /// the vanilla partial_reduce function. Specifically it implements `n = n0 + (2^32 -1)n1`, which
    /// is equivalent to `n = n0 + (2^32 -1)n1 - n2` if and only if n2 == 0
    ///
    /// The correctness of this functions follows from this equation in section 2 of
    /// https://github.com/mir-protocol/plonky2/blob/main/plonky2.pdf
    pub fn reduce_u96(low_order_bits: u64, high_order_bits: u32) -> u64 {
        let n0 = low_order_bits;
        let n1: u64 = ((high_order_bits as u64) << 32) - (high_order_bits as u64);
        let (mut result, overflow) = n0.overflowing_add(n1);
        if overflow {
            // If we overflowed, this subtraction will underflow yielding the correct result
            result = Self::sub_p_from(result)
        }
        result
    }

    /// Reduce the result of u64::carrying_add to a u64.
    /// The correctness of this reduction follows from the equation `n = n0 + (2^32 -1)n1 - n2`.
    /// Because addition can never overflow by more than a single bit, n2 - 0. And n1 is either 1
    /// (if the addition overflowed) or 0 (if it did not overflow).
    /// Now we simply hard code the results of each case:
    /// `n = n0 + (2^32 -1)` if overflowed
    /// else: `n = n0`
    pub fn reduce_addition(overflowing_add_result: (u64, bool)) -> u64 {
        // let (result, overflowed) = overflowing_add_result;
        // match overflowed {
        //     // If we overflowed, compute n0 + (2^32 -1) = n0 + U32_MAX
        //     true => {
        //         // If the result was already greater than P, subtract P before adding u32_MAX
        //         if result >= P {
        //             // Equivalent to sub_p_from(result) + U32_MAX.
        //             // sub_p_from(result) = lhs.wrapping_add(U32_MAX). Thus the wrapped portion of
        //             // result is less than P, so we can safely add U32_MAX without overflowing again.
        //             // For speed, we simply combine the two additions.
        //             return result.wrapping_add(U32_MAX_DOUBLED);
        //         }
        //         // Otherwise, result < P => result < (2^64 - 2^32 +1)
        //         // => result + (2^32 - 1) < (2^64 - 2^32 +1) (2^32 - 1)
        //         // => result + (2^32 - 1)  < 2^64. Therefore, the operation cannot overflow
        //         return unsafe { result.unchecked_add(U32_MAX) };
        //     }
        //     false => result,
        // }
        if overflowing_add_result.1 {
            return Self::reduce_u96(overflowing_add_result.0, 1);
        }
        overflowing_add_result.0
    }

    #[inline]
    // Taken from section 2 of https://github.com/mir-protocol/plonky2/blob/main/plonky2.pdf
    fn sub_p_from(lhs: u64) -> u64 {
        lhs.wrapping_add(U32_MAX)
    }
}

impl Uint for u64 {
    fn is_zero(&self) -> bool {
        *self == 0
    }

    fn zero() -> Self {
        0
    }

    fn one() -> Self {
        1
    }

    fn is_even(&self) -> bool {
        *self & 1 == 0
    }

    fn is_odd(&self) -> bool {
        *self & 1 == 1
    }
    fn to_bytes_le(&self, out: &mut [u8]) {
        out.clone_from_slice(&self.to_le_bytes())
    }
}

impl PrimeField for GoldilocksField {
    type Elem = u64;

    fn reduce(n: Self::Elem) -> Self::Elem {
        if n >= P {
            return n - P;
        }
        n
    }

    #[inline]
    fn p() -> Self::Elem {
        P
    }

    fn generator() -> Self::Elem {
        unimplemented!()
    }

    fn primitive_nth_root(_n: Self::Elem) -> Result<Self::Elem, crate::field::NoNthRootError> {
        unimplemented!()
    }

    fn sample(random: &[u8]) -> Self::Elem {
        let random: Vec<u8> = if random.len() < 8 {
            random
                .iter()
                .map(|v| *v)
                .chain(std::iter::repeat(0))
                .take(8)
                .collect()
        } else {
            random.iter().map(|v| *v).take(8).collect()
        };
        let bytes = random.try_into().unwrap();
        u64::from_be_bytes(bytes)
    }

    fn zero() -> Self::Elem {
        Uint::zero()
    }

    fn one() -> Self::Elem {
        Uint::one()
    }

    fn add(x: u64, y: u64) -> u64 {
        Self::reduce_addition(x.overflowing_add(y))
    }

    fn subtract(lhs: u64, rhs: u64) -> u64 {
        let rhs = if rhs >= P {
            GoldilocksField::sub_p_from(rhs)
        } else {
            rhs
        };
        let (result, underflow) = lhs.overflowing_sub(rhs);
        if underflow {
            // If we overflowed, this subtraction will underflow yielding the correct result
            return result.wrapping_add(P);
        }
        result
    }

    fn multiply(lhs: u64, rhs: u64) -> u64 {
        Self::reduce_multiplication(lhs.widening_mul(rhs))
    }

    fn negate(value: Self::Elem) -> Self::Elem {
        Self::p() - Self::reduce(value)
    }

    fn xgcd(lhs: u64, rhs: u64) -> (u64, u64, u64) {
        let (mut old_r, mut r) = (lhs, rhs);
        let (mut old_s, mut s) = (1, 0);
        let (mut old_t, mut t) = (0, 1);

        // if r >= P {
        //     r = GoldilocksField::sub_p_from(r)
        // }
        while r != 0 {
            let quotient = old_r / r;
            (old_r, r) = (r, Self::subtract(old_r, Self::multiply(quotient, r)));
            (old_s, s) = (s, Self::subtract(old_s, Self::multiply(quotient, s)));
            (old_t, t) = (t, Self::subtract(old_t, Self::multiply(quotient, t)));
            if r >= P {
                r = GoldilocksField::sub_p_from(r)
            }
        }
        (old_s, old_t, old_r)
    }

    fn inverse(value: Self::Elem) -> Self::Elem {
        let (a, _, _) = Self::xgcd(value, Self::p());
        a
    }

    fn divide(l: Self::Elem, r: Self::Elem) -> Result<Self::Elem, crate::field::DivByZeroError> {
        if r.is_zero() {
            return Err(crate::field::DivByZeroError);
        }

        Ok(Self::multiply(l, Self::inverse(r)))
    }
}

// pub trait GoldilocksMath:
//     Clone
//     + std::fmt::Debug
//     + Copy
//     + Add<Self>
//     + Mul<Self>
//     + Div<Self>
//     + Sub<Self>
//     + PartialOrd
//     + PartialEq
// {
//     fn inverse(self) -> Self;
//     fn xgcd(self, rhs: Self) -> (Self, Self, Self);

//     fn pow(self, rhs: u64) -> Self;
//     // fn zero() -> Self;
//     // fn one() -> Self;
//     // fn zeros(len: usize) -> Vec<Self>;
// }

#[cfg(test)]
mod tests {
    use crate::{field::PrimeField, GoldilocksField, P};

    #[test]
    fn test_add() {
        assert_eq!(GoldilocksField::add(P, P), P);
        assert_eq!(GoldilocksField::add(P, P / 2), P / 2);
        assert_eq!(GoldilocksField::add(P + 1, P + 1), P + 2);
        assert_eq!(GoldilocksField::add(u64::MAX, u64::MAX), 8589934588)
    }

    #[test]
    fn test_sub() {
        assert_eq!(GoldilocksField::subtract(P, P), P);
        assert_eq!(GoldilocksField::subtract(P, P / 2), (P / 2) + 1);
        assert_eq!(GoldilocksField::subtract(P + 1, P + 1), P);
        assert_eq!(GoldilocksField::subtract(u64::MAX, u64::MAX), P);
        assert_eq!(GoldilocksField::subtract(0, u64::MAX), 18446744065119617027);
        assert_eq!(GoldilocksField::subtract(0, P - 2), 2);
    }
}

#[cfg(test)]
mod integ_tests {
    use crate::{
        field::PrimeField, field_elem::FieldElement, uint::Uint, GoldilocksElement,
        P as GoldilocksPrime,
    };
    use proptest::prelude::*;

    #[derive(Debug, Clone, Copy)]
    struct DefaultGoldilocksField;

    impl Uint for u128 {
        fn is_zero(&self) -> bool {
            *self == 0
        }

        fn zero() -> Self {
            0
        }

        fn one() -> Self {
            1
        }

        fn is_even(&self) -> bool {
            *self & 1 == 0
        }

        fn is_odd(&self) -> bool {
            *self & 1 == 0
        }

        fn to_bytes_le(&self, out: &mut [u8]) {
            out.clone_from_slice(&self.to_le_bytes())
        }
    }

    impl PrimeField for DefaultGoldilocksField {
        type Elem = u128;

        fn reduce(n: Self::Elem) -> Self::Elem {
            n % Self::p()
        }

        fn p() -> Self::Elem {
            GoldilocksPrime as u128
        }

        fn generator() -> Self::Elem {
            todo!()
        }

        fn primitive_nth_root(_n: Self::Elem) -> Result<Self::Elem, crate::field::NoNthRootError> {
            todo!()
        }

        fn sample(_random: &[u8]) -> Self::Elem {
            todo!()
        }
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(100000))]
        #[test]
        fn adds_match(a: u64, b: u64) {
            let goldilocks_result = GoldilocksElement::from_u64(a) + GoldilocksElement::from_u64(b);
            let default_result = DefaultGoldilocksField::add(a as u128, b as u128);
            prop_assert_eq!(goldilocks_result, default_result as u64)
        }

        #[test]
        fn subtracts_match(a: u64, b: u64) {
            let goldilocks_result = GoldilocksElement::from_u64(a) - GoldilocksElement::from_u64(b);
            let default_result = DefaultGoldilocksField::subtract(a as u128, b as u128);
            prop_assert_eq!(goldilocks_result, default_result as u64)
        }


        #[test]
        fn multiplies_match(a: u64, b: u64) {
            let goldilocks_result = GoldilocksElement::from_u64(a) * GoldilocksElement::from_u64(b);
            let default_result = DefaultGoldilocksField::multiply(a as u128, b as u128);
            prop_assert_eq!(goldilocks_result, default_result as u64)
        }

        #[test]
        fn inverts_match(a: u64) {
            let goldilocks_inverse = GoldilocksElement::from_u64(a).inverse();

            let default_result = DefaultGoldilocksField::inverse(a as u128);
            prop_assert_eq!(goldilocks_inverse,  default_result as u64)
        }

        #[test]
        fn divides_match(a: u64, b: u64) {
            let div_result = GoldilocksElement::from_u64(a) / GoldilocksElement::from_u64(b);
            let default_result = DefaultGoldilocksField::divide(a as u128, b as u128);
            prop_assert_eq!(div_result, default_result.unwrap_or(0) as u64)
        }

        #[test]
        fn negates_match(a: u64) {
            let neg_result = -GoldilocksElement::from_u64(a);
            let default_result = DefaultGoldilocksField::negate(a as u128);
            prop_assert_eq!(neg_result, default_result as u64)
        }

        #[test]
        fn pows_match(a: u64, b: usize) {
            let pow_result  = GoldilocksElement::from_u64(a).pow(b);
            let default_result = DefaultGoldilocksField::pow(a as u128, b);
            prop_assert_eq!(pow_result, default_result as u64)
        }
    }
}
