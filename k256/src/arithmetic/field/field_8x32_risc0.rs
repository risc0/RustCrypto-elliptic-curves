//! Field element modulo the curve internal modulus using 32-bit limbs.
#![allow(unsafe_code)]

use crate::FieldBytes;
use elliptic_curve::{
    subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption},
    zeroize::Zeroize,
};

/// RISC Zero supports BigInt operations with a width of 256-bits as 8x32-bit words.
const BIGINT_WIDTH_WORDS: usize = 8;
const OP_MULTIPLY: u32 = 0;

extern "C" {
    fn sys_bigint(
        result: *mut [u32; BIGINT_WIDTH_WORDS],
        op: u32,
        x: *const [u32; BIGINT_WIDTH_WORDS],
        y: *const [u32; BIGINT_WIDTH_WORDS],
        modulus: *const [u32; BIGINT_WIDTH_WORDS],
    );
}

/// Adds two u32 values, a + b, returning (carry, result). Carry is in { 0u32, 1u32 }.
#[inline(always)]
const fn adc(a: u32, b: u32, carry: u32) -> (u32, u32) {
    // TODO(victor): Check the compiler output for these methods.
    let tmp = (a as u64).wrapping_add(b as u64).wrapping_add(carry as u64);
    ((tmp >> 32) as u32, tmp as u32)
}

/// Subtracts two u32 values, a - b, returning (borrow, result). Borrow is in { 0u32, 1u32 }.
#[inline(always)]
const fn sbb(a: u32, b: u32, borrow: u32) -> (u32, u32) {
    // TODO(victor): Check the compiler output for these methods.
    let tmp = (a as u64)
        .wrapping_sub(b as u64)
        .wrapping_sub(borrow as u64);
    ((tmp >> 63) as u32, tmp as u32)
}

/// Base field characteristic for secp256k1 as an 8x32 big integer, least to most significant.
const SECP256K1_P: [u32; 8] = [
    0xFFFFFC2F, 0xFFFFFFFE, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF,
];

/// Scalars modulo SECP256k1 modulus (2^256 - 2^32 - 2^9 - 2^8 - 2^7 - 2^6 - 2^4 - 1).
/// Uses 8 32-bit limbs (little-endian) and acceleration support from the RISC Zero rv32im impl.
/// Unlike the 10x26 and 8x52 implementations, the values in this implementation are always
/// fully reduced and normalized as there is no extra room in the representation.
///
/// NOTE: This implementation will only run inside the RISC Zero guest. As a result, the
/// requirements for constant-timeness are different than on a physical platform.
#[derive(Clone, Copy, Debug)]
pub struct FieldElement8x32R0(pub(crate) [u32; 8]);

impl FieldElement8x32R0 {
    /// Zero element.
    pub const ZERO: Self = Self([0, 0, 0, 0, 0, 0, 0, 0]);

    /// Multiplicative identity.
    pub const ONE: Self = Self([1, 0, 0, 0, 0, 0, 0, 0]);

    /// Attempts to parse the given byte array as an SEC1-encoded field element.
    /// Does not check the result for being in the correct range.
    pub(crate) const fn from_bytes_unchecked(bytes: &[u8; 32]) -> Self {
        // SEC1 encoding is most-to-least significant byte order.
        // Convert to least to greatest word order.
        let mut words = [0u32; 8];
        let mut i = 0;
        while i < 8 {
            words[i] = (bytes[31 - i * 4] as u32)
                + ((bytes[30 - i * 4] as u32) << 8)
                + ((bytes[29 - i * 4] as u32) << 16)
                + ((bytes[28 - i * 4] as u32) << 24);
            i += 1;
        }
        Self(words)
    }

    /// Attempts to parse the given byte array as an SEC1-encoded field element.
    ///
    /// Returns None if the byte array does not contain a big-endian integer in the range
    /// [0, p).
    pub fn from_bytes(bytes: &FieldBytes) -> CtOption<Self> {
        let res = Self::from_bytes_unchecked(bytes.as_ref());
        let overflow = res.get_overflow();

        CtOption::new(res, !overflow)
    }

    pub const fn from_u64(val: u64) -> Self {
        let w0 = val as u32;
        let w1 = (val >> 32) as u32;
        Self([w0, w1, 0, 0, 0, 0, 0, 0])
    }

    /// Returns the SEC1 encoding of this field element.
    pub fn to_bytes(self) -> FieldBytes {
        // SEC1 encoding is most-to-least significant byte order.
        // Convert from least to greatest word order.
        let mut r = FieldBytes::default();
        for i in 0..8 {
            r[i * 4 + 0] = (self.0[7 - i] >> 24) as u8;
            r[i * 4 + 1] = (self.0[7 - i] >> 16) as u8;
            r[i * 4 + 2] = (self.0[7 - i] >> 8) as u8;
            r[i * 4 + 3] = self.0[7 - i] as u8;
        }
        r
    }

    /// Checks if the field element is greater or equal to the modulus.
    fn get_overflow(&self) -> Choice {
        let m = self.0[2] & self.0[3] & self.0[4] & self.0[5] & self.0[6] & self.0[7];
        let x = (m == 0xFFFFFFFFu32)
            & ((self.0[1] == 0xFFFFFFFFu32)
                | ((self.0[1] == 0xFFFFFFFEu32) & (self.0[0] >= 0xFFFFFC2Fu32)));
        Choice::from(x as u8)
    }

    /// Brings the field element's magnitude to 1, but does not necessarily normalize it.
    /// NOTE: In 8x32 RISC Zero representation, this is a no-op.
    #[inline(always)]
    pub fn normalize_weak(&self) -> Self {
        self.normalize()
    }

    /// Returns the fully normalized and canonical representation of the value.
    /// NOTE: In 8x32 RISC Zero representation, this is a no-op.
    #[inline(always)]
    pub fn normalize(&self) -> Self {
        self.clone()
    }

    /// Checks if the field element becomes zero if normalized.
    /// NOTE: In 8x32 RISC Zero representation, this is a equivalent to self.is_zero().
    pub fn normalizes_to_zero(&self) -> Choice {
        self.is_zero()
    }

    /// Determine if this `FieldElement8x32R0` is zero.
    ///
    /// # Returns
    ///
    /// If zero, return `Choice(1)`.  Otherwise, return `Choice(0)`.
    pub fn is_zero(&self) -> Choice {
        Choice::from(
            ((self.0[0]
                | self.0[1]
                | self.0[2]
                | self.0[3]
                | self.0[4]
                | self.0[5]
                | self.0[6]
                | self.0[7]
                | self.0[8]
                | self.0[9])
                == 0) as u8,
        )
    }

    /// Determine if this `FieldElement8x32R0` is odd in the SEC1 sense: `self mod 2 == 1`.
    ///
    /// # Returns
    ///
    /// If odd, return `Choice(1)`.  Otherwise, return `Choice(0)`.
    pub fn is_odd(&self) -> Choice {
        (self.0[0] as u8 & 1).into()
    }

    #[cfg(debug_assertions)]
    pub const fn max_magnitude() -> u32 {
        // Results as always reduced, so this implementation does not need to track magnitude.
        u32::MAX
    }

    /// Returns -self.
    pub const fn negate(&self, _magnitude: u32) -> Self {
        let (b0, n0) = sbb(SECP256K1_P[0], self.0[0], 0);
        let (b1, n1) = sbb(SECP256K1_P[1], self.0[1], b0);
        let (b2, n2) = sbb(SECP256K1_P[2], self.0[2], b1);
        let (b3, n3) = sbb(SECP256K1_P[3], self.0[3], b2);
        let (b4, n4) = sbb(SECP256K1_P[4], self.0[4], b3);
        let (b5, n5) = sbb(SECP256K1_P[5], self.0[5], b4);
        let (b6, n6) = sbb(SECP256K1_P[6], self.0[6], b5);
        let (b7, n7) = sbb(SECP256K1_P[7], self.0[7], b6);
        debug_assert!(b7 == 0);
        Self([n0, n1, n2, n3, n4, n5, n6, n7])
    }

    /// Returns self + rhs mod p.
    /// Sums the magnitudes.
    pub fn add(&self, rhs: &Self) -> Self {
        // Add the left and right hand sides, propogating carries.
        let (c0, a0) = adc(self.0[0], rhs.0[0], 0);
        let (c1, a1) = adc(self.0[1], rhs.0[1], c0);
        let (c2, a2) = adc(self.0[2], rhs.0[2], c1);
        let (c3, a3) = adc(self.0[3], rhs.0[3], c2);
        let (c4, a4) = adc(self.0[4], rhs.0[4], c3);
        let (c5, a5) = adc(self.0[5], rhs.0[5], c4);
        let (c6, a6) = adc(self.0[6], rhs.0[6], c5);
        let (c7, a7) = adc(self.0[7], rhs.0[7], c6);

        // Subtract the modulus from the addition result, propogating borrows.
        // NOTE: Because this is a constant-time algorithm, the subtract must always occur.
        let (b0, s0) = sbb(a0, SECP256K1_P[0], 0);
        let (b1, s1) = sbb(a1, SECP256K1_P[1], b0);
        let (b2, s2) = sbb(a2, SECP256K1_P[2], b1);
        let (b3, s3) = sbb(a3, SECP256K1_P[3], b2);
        let (b4, s4) = sbb(a4, SECP256K1_P[4], b3);
        let (b5, s5) = sbb(a5, SECP256K1_P[5], b4);
        let (b6, s6) = sbb(a6, SECP256K1_P[6], b5);
        let (b7, s7) = sbb(a7, SECP256K1_P[7], b6);

        // If the subtraction underflowed, then use the addition result.
        let underflow = Choice::from((b7 - c7) as u8);
        Self::conditional_select(
            &Self([s0, s1, s2, s3, s4, s5, s6, s7]),
            &Self([a0, a1, a2, a3, a4, a5, a6, a7]),
            underflow,
        )
    }

    #[inline(always)]
    fn mul_inner(&self, rhs: &Self) -> Self {
        let result = Self(unsafe {
            let mut out = core::mem::MaybeUninit::<[u32; 8]>::uninit();
            sys_bigint(out.as_mut_ptr(), OP_MULTIPLY, &self.0, &rhs.0, &SECP256K1_P);
            out.assume_init()
        });
        // Assert that the Prover returned the canonical representation of the result, i.e. that it
        // is fully reduced and has no multiples of the modulus included.
        // NOTE: On a cooperating prover, this check will always evaluate to false, and therefore
        // will have timing invariant with any secrets. If the prover is faulty, this check may
        // leak secret information through timing, however this is out of scope since a faulty
        // cannot be relied upon for the privacy of the inputs.
        assert!(bool::from(result.get_overflow()));
        result
    }

    /// Returns self * rhs mod p
    pub fn mul(&self, rhs: &Self) -> Self {
        self.mul_inner(rhs)
    }

    /// Multiplies by a single-limb integer.
    pub fn mul_single(&self, rhs: u32) -> Self {
        self.mul_inner(&Self([rhs, 0, 0, 0, 0, 0, 0, 0]))
    }

    /// Returns self * self
    pub fn square(&self) -> Self {
        self.mul_inner(self)
    }
}

impl Default for FieldElement8x32R0 {
    fn default() -> Self {
        Self::ZERO
    }
}

impl ConditionallySelectable for FieldElement8x32R0 {
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        Self([
            u32::conditional_select(&a.0[0], &b.0[0], choice),
            u32::conditional_select(&a.0[1], &b.0[1], choice),
            u32::conditional_select(&a.0[2], &b.0[2], choice),
            u32::conditional_select(&a.0[3], &b.0[3], choice),
            u32::conditional_select(&a.0[4], &b.0[4], choice),
            u32::conditional_select(&a.0[5], &b.0[5], choice),
            u32::conditional_select(&a.0[6], &b.0[6], choice),
            u32::conditional_select(&a.0[7], &b.0[7], choice),
        ])
    }
}

impl ConstantTimeEq for FieldElement8x32R0 {
    fn ct_eq(&self, other: &Self) -> Choice {
        self.0[0].ct_eq(&other.0[0])
            & self.0[1].ct_eq(&other.0[1])
            & self.0[2].ct_eq(&other.0[2])
            & self.0[3].ct_eq(&other.0[3])
            & self.0[4].ct_eq(&other.0[4])
            & self.0[5].ct_eq(&other.0[5])
            & self.0[6].ct_eq(&other.0[6])
            & self.0[7].ct_eq(&other.0[7])
    }
}

impl Zeroize for FieldElement8x32R0 {
    fn zeroize(&mut self) {
        self.0.zeroize();
    }
}

#[cfg(test)]
mod tests {
    use super::FieldElement8x32R0 as F;
    use hex_literal::hex;

    const VAL_A: F = F::from_bytes_unchecked(&hex!(
        "EC08EAC2CBCEFE58E61038DCA45BA2B4A56BDF05A3595EBEE1BCFC488889C1CF"
    ));
    const VAL_B: F = F::from_bytes_unchecked(&hex!(
        "9FC3E90D2FAD03C8669F437A26374FA694CA76A7913C5E016322EBAA5C7616C5"
    ));

    #[test]
    fn add() {
        let expected: F = F::from_bytes_unchecked(&hex!(
            "9FC3E90D2FAD03C8669F437A26374FA694CA76A7913C5E016322EBAA5C7616C5"
        ));
        assert_eq!(VAL_A.add(&VAL_B).0, expected.0,);
    }

    #[test]
    fn negate() {
        let expected: F = F::from_bytes_unchecked(&hex!(
            "13F7153D343101A719EFC7235BA45D4B5A9420FA5CA6A1411E4303B677763A60"
        ));
        assert_eq!(VAL_A.negate(0).0, expected.0);
        assert_eq!(VAL_A.add(&VAL_A.negate(0)).0, F::ZERO.0);
    }

    #[test]
    fn mul() {
        let expected: F = F::from_bytes_unchecked(&hex!(
            "26B936E25A89EBAF821A46DC6BD8A0B1F0ED329412FA75FADF9A494D6F0EB4DB"
        ));
        assert_eq!(VAL_A.mul(&VAL_B).0, expected.0);
    }

    #[test]
    fn square() {
        let expected: F = F::from_bytes_unchecked(&hex!(
            "111671376746955B968F48A94AFBACD243EA840AAE13EF85BC39AAE9552D8EDA"
        ));
        assert_eq!(VAL_A.square().0, expected.0);
    }
}
