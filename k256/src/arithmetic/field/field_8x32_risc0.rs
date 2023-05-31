//! Field element modulo the curve internal modulus using 32-bit limbs.
#![allow(unsafe_code)]

use crate::FieldBytes;
use elliptic_curve::{
    bigint::{risc0, ArrayEncoding, Integer, Limb, Zero, U256},
    subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption},
    zeroize::Zeroize,
};

/// Base field characteristic for secp256k1 as an 8x32 big integer, least to most significant.
const MODULUS: U256 =
    U256::from_be_hex("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F");

/// Low two words of 2^256 - MODULUS, used for correcting the value after addition mod 2^256.
const MODULUS_CORRECTION_0: u32 = U256::ZERO.wrapping_sub(&MODULUS).as_words()[0];
const MODULUS_CORRECTION_1: u32 = U256::ZERO.wrapping_sub(&MODULUS).as_words()[1];

/// Scalars modulo SECP256k1 modulus (2^256 - 2^32 - 2^9 - 2^8 - 2^7 - 2^6 - 2^4 - 1).
/// Uses 8 32-bit limbs (little-endian) and acceleration support from the RISC Zero rv32im impl.
/// Unlike the 10x26 and 8x52 implementations, the values in this implementation are always
/// fully reduced and normalized as there is no extra room in the representation.
///
/// NOTE: This implementation will only run inside the RISC Zero guest. As a result, the
/// requirements for constant-timeness are different than on a physical platform.
#[derive(Clone, Copy, Debug)]
pub struct FieldElement8x32R0(pub(crate) U256);

impl FieldElement8x32R0 {
    /// Zero element.
    pub const ZERO: Self = Self(U256::ZERO);

    /// Multiplicative identity.
    pub const ONE: Self = Self(U256::ONE);

    /// Attempts to parse the given byte array as an SEC1-encoded field element.
    /// Does not check the result for being in the correct range.
    pub(crate) const fn from_bytes_unchecked(bytes: &[u8; 32]) -> Self {
        Self(U256::from_be_slice(&bytes.as_slice()))
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
        Self(U256::from_words([w0, w1, 0, 0, 0, 0, 0, 0]))
    }

    pub const fn from_i64(val: i64) -> Self {
        // Compute val_abs = |val|
        let val_mask = val >> 63;
        let val_abs = ((val + val_mask) ^ val_mask) as u64;

        Self::from_u64(val_abs).negate_const()
    }

    /// Returns the SEC1 encoding of this field element.
    pub fn to_bytes(self) -> FieldBytes {
        self.0.to_be_byte_array()
    }

    /// Checks if the field element is greater or equal to the modulus.
    fn get_overflow(&self) -> Choice {
        let words = self.0.as_words();
        let m = words[2] & words[3] & words[4] & words[5] & words[6] & words[7];
        let x = (m == 0xFFFFFFFFu32)
            & ((words[1] == 0xFFFFFFFFu32)
                | ((words[1] == 0xFFFFFFFEu32) & (words[0] >= 0xFFFFFC2Fu32)));
        Choice::from(x as u8)
    }

    /// Brings the field element's magnitude to 1, but does not necessarily normalize it.
    ///
    /// NOTE: In RISC Zero, this is a no-op since weak normalization is not an operation that
    /// needs to be performed between calls to arithmetic routines.
    #[inline(always)]
    pub fn normalize_weak(&self) -> Self {
        self.clone()
    }

    /// Returns the fully normalized and canonical representation of the value.
    #[inline(always)]
    pub fn normalize(&self) -> Self {
        // When the prover is cooperative, the value is always normalized.
        assert!(!bool::from(self.get_overflow()));
        self.clone()
    }

    /// Checks if the field element becomes zero if normalized.
    pub fn normalizes_to_zero(&self) -> Choice {
        self.0.ct_eq(&U256::ZERO) | self.0.ct_eq(&MODULUS)
    }

    /// Determine if this `FieldElement8x32R0` is zero.
    ///
    /// # Returns
    ///
    /// If zero, return `Choice(1)`.  Otherwise, return `Choice(0)`.
    pub fn is_zero(&self) -> Choice {
        self.0.is_zero()
    }

    /// Determine if this `FieldElement8x32R0` is odd in the SEC1 sense: `self mod 2 == 1`.
    ///
    /// Value must be normalized before calling is_odd.
    ///
    /// # Returns
    ///
    /// If odd, return `Choice(1)`.  Otherwise, return `Choice(0)`.
    pub fn is_odd(&self) -> Choice {
        self.0.is_odd()
    }

    #[cfg(debug_assertions)]
    pub const fn max_magnitude() -> u32 {
        // Results as always reduced, so this implementation does not need to track magnitude.
        u32::MAX
    }

    /// Returns -self.
    const fn negate_const(&self) -> Self {
        let (s, borrow) = MODULUS.sbb(&self.0, Limb(0));
        assert!(borrow.0 == 0);
        Self(s)
    }

    /// Returns -self.
    pub fn negate(&self, _magnitude: u32) -> Self {
        self.mul(&Self::ONE.negate_const())
    }

    /// Returns self + rhs mod p.
    /// Sums the magnitudes.
    pub fn add(&self, rhs: &Self) -> Self {
        let (a, carry) = self.0.adc(&rhs.0, Limb(0));

        // If a carry or overflow of the modulus occured, we need to add 2^256 - p.
        // c0 and c1 and the two non-zero limbs of the correction value.
        let denorm = Self(a).get_overflow().unwrap_u8() as u32;
        let mask = carry.0 | denorm;
        let c0 = MODULUS_CORRECTION_0 * mask;
        let c1 = MODULUS_CORRECTION_1 * mask;
        let correction = U256::from_words([c0, c1, 0, 0, 0, 0, 0, 0]);

        Self(a.wrapping_add(&correction))
    }

    /// Returns self * rhs mod p
    pub fn mul(&self, rhs: &Self) -> Self {
        Self(risc0::modmul_u256_denormalized(&self.0, &rhs.0, &MODULUS))
    }

    /// Multiplies by a single-limb integer.
    pub fn mul_single(&self, rhs: u32) -> Self {
        Self(risc0::modmul_u256_denormalized(
            &self.0,
            &U256::from_words([rhs, 0, 0, 0, 0, 0, 0, 0]),
            &MODULUS,
        ))
    }

    /// Returns self * self
    pub fn square(&self) -> Self {
        Self(risc0::modmul_u256_denormalized(&self.0, &self.0, &MODULUS))
    }
}

impl Default for FieldElement8x32R0 {
    fn default() -> Self {
        Self::ZERO
    }
}

impl ConditionallySelectable for FieldElement8x32R0 {
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        Self(U256::conditional_select(&a.0, &b.0, choice))
    }
}

impl ConstantTimeEq for FieldElement8x32R0 {
    fn ct_eq(&self, other: &Self) -> Choice {
        self.0.ct_eq(&other.0)
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

    extern crate alloc;

    fn as_hex(&elem: &F) -> alloc::string::String {
        ::hex::encode_upper(elem.to_bytes())
    }

    #[test]
    fn add() {
        let expected: F = F::from_bytes_unchecked(&hex!(
            "8BCCD3CFFB7C02214CAF7C56CA92F25B3A3655AD3495BCC044DFE7F3E4FFDC65"
        ));
        assert_eq!(as_hex(&VAL_A.add(&VAL_B)), as_hex(&expected));
    }

    // Tests the other "code path" returning the reduced or non-reduced result.
    #[test]
    fn add_negated() {
        let expected: F = F::from_bytes_unchecked(&hex!(
            "74332C300483FDDEB35083A9356D0DA4C5C9AA52CB6A433FBB20180B1B001FCA"
        ));
        assert_eq!(
            as_hex(&VAL_A.negate(0).add(&VAL_B.negate(0))),
            as_hex(&expected)
        );
    }

    #[test]
    fn negate() {
        let expected: F = F::from_bytes_unchecked(&hex!(
            "13F7153D343101A719EFC7235BA45D4B5A9420FA5CA6A1411E4303B677763A60"
        ));
        assert_eq!(as_hex(&VAL_A.negate(0)), as_hex(&expected));
        assert_eq!(as_hex(&VAL_A.add(&VAL_A.negate(0))), as_hex(&F::ZERO));
    }

    #[test]
    fn mul() {
        let expected: F = F::from_bytes_unchecked(&hex!(
            "26B936E25A89EBAF821A46DC6BD8A0B1F0ED329412FA75FADF9A494D6F0EB4DB"
        ));
        assert_eq!(as_hex(&VAL_A.mul(&VAL_B)), as_hex(&expected));
    }

    #[test]
    fn mul_zero() {
        assert_eq!(as_hex(&VAL_A.mul(&F::ZERO)), as_hex(&F::ZERO));
        assert_eq!(as_hex(&VAL_B.mul(&F::ZERO)), as_hex(&F::ZERO));
        assert_eq!(as_hex(&F::ZERO.mul(&F::ZERO)), as_hex(&F::ZERO));
        assert_eq!(as_hex(&F::ONE.mul(&F::ZERO)), as_hex(&F::ZERO));
        assert_eq!(as_hex(&F::ONE.negate(0).mul(&F::ZERO)), as_hex(&F::ZERO));
    }

    #[test]
    fn mul_one() {
        assert_eq!(as_hex(&VAL_A.mul(&F::ONE)), as_hex(&VAL_A));
        assert_eq!(as_hex(&VAL_B.mul(&F::ONE)), as_hex(&VAL_B));
        assert_eq!(as_hex(&F::ZERO.mul(&F::ONE)), as_hex(&F::ZERO));
        assert_eq!(as_hex(&F::ONE.mul(&F::ONE)), as_hex(&F::ONE));
        assert_eq!(
            as_hex(&F::ONE.negate(0).mul(&F::ONE)),
            as_hex(&F::ONE.negate(0))
        );
    }

    #[test]
    fn square() {
        let expected: F = F::from_bytes_unchecked(&hex!(
            "111671376746955B968F48A94AFBACD243EA840AAE13EF85BC39AAE9552D8EDA"
        ));
        assert_eq!(as_hex(&VAL_A.square()), as_hex(&expected));
    }
}
