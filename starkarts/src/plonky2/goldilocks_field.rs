use std::ops::{Add, Div, Mul, Sub};

#[derive(Debug, Clone, Copy)]
pub struct GoldilocksField;
pub const P: u64 = u64::MAX - (u32::MAX as u64) + 1;
const LOW_32_BITS: u64 = std::u32::MAX as u64;
const U32_MAX: u64 = std::u32::MAX as u64;
const U32_MAX_DOUBLED: u64 = (std::u32::MAX as u64) * 2;

///  A 64 bit field allowing very arithmetic on modern CPUs.
///  https://github.com/mir-protocol/plonky2/blob/main/plonky2.pdf
impl GoldilocksField {
    /// A (somewhat) optimized function to compute X % p for any 128 bit integer x,
    /// Based on section 2 of https://github.com/mir-protocol/plonky2/blob/main/plonky2.pdf
    /// Note that this function computes a reduction to a u64, but does not put the result into
    /// canonical form. Specifically, the result may be larger than P  so a single subtraction
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
        let n1: u64 = (high_order_bits as u64) << 32 - high_order_bits;
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
        let (result, overflowed) = overflowing_add_result;
        match overflowed {
            // If we overflowed, compute n0 + (2^32 -1) = n0 + U32_MAX
            true => {
                // If the result was already greater than P, subtract P before adding u32_MAX
                if result >= P {
                    // Equivalent to sub_p_from(result) + U32_MAX.
                    // sub_p_from(result) = lhs.wrapping_add(U32_MAX). Thus the wrapped portion of
                    // result is less than P, so we can safely add U32_MAX without overflowing again.
                    // For speed, we simply combine the two additions.
                    return result.wrapping_add(U32_MAX_DOUBLED);
                }
                // Otherwise, result < P => result < (2^64 - 2^32 +1)
                // => result + (2^32 - 1) < (2^64 - 2^32 +1) (2^32 - 1)
                // => result + (2^32 - 1)  < 2^64. Therefore, the operation cannot overflow
                return unsafe { result.unchecked_add(U32_MAX) };
            }
            false => result,
        }
    }

    pub fn add(x: u64, y: u64) -> u64 {
        Self::reduce_addition(x.overflowing_add(y))
    }

    pub fn xgcd(lhs: u64, rhs: u64) -> (u64, u64, u64) {
        let (mut old_r, mut r) = (lhs, rhs);
        let (mut old_s, mut s) = (1, 0);
        let (mut old_t, mut t) = (0, 1);

        while r != 0 {
            let quotient = old_r / r;
            (old_r, r) = (r, Self::subtract(old_r, Self::multiply(quotient, r)));
            (old_s, s) = (s, Self::subtract(old_s, Self::multiply(quotient, s)));
            (old_t, t) = (t, Self::subtract(old_t, Self::multiply(quotient, t)));
        }
        (old_s, old_t, old_r)
    }

    pub fn subtract(lhs: u64, rhs: u64) -> u64 {
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

    pub fn multiply(lhs: u64, rhs: u64) -> u64 {
        Self::reduce_multiplication(lhs.widening_mul(rhs))
    }

    #[inline]
    // Taken from section 2 of https://github.com/mir-protocol/plonky2/blob/main/plonky2.pdf
    fn sub_p_from(lhs: u64) -> u64 {
        lhs.wrapping_add(U32_MAX)
    }
}

pub trait GoldilocksMath:
    Clone
    + std::fmt::Debug
    + Copy
    + Add<Self>
    + Mul<Self>
    + Div<Self>
    + Sub<Self>
    + PartialOrd
    + PartialEq
{
    fn inverse(self) -> Self;
    fn xgcd(self, rhs: Self) -> (Self, Self, Self);

    fn pow(self, rhs: u64) -> Self;
    // fn zero() -> Self;
    // fn one() -> Self;
    // fn zeros(len: usize) -> Vec<Self>;
}

#[cfg(test)]
mod tests {
    use crate::{GoldilocksField, P};

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