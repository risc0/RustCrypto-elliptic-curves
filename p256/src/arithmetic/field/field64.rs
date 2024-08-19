//! 64-bit secp256r1 field element algorithms.

use super::{FieldBytes, FieldElement, MODULUS};
use crate::arithmetic::util::{adc, mac, sbb, u256_to_u64x4, u64x4_to_u256};
use elliptic_curve::{
    bigint::{ArrayEncoding, U256},
    rand_core::RngCore,
    subtle::{Choice, ConstantTimeEq},
};

/// R = 2^256 mod p
const R: FieldElement = FieldElement(U256::from_be_hex(
    "00000000fffffffeffffffffffffffffffffffff000000000000000000000001",
));

/// R^2 = 2^512 mod p
const R2: FieldElement = FieldElement(U256::from_be_hex(
    "00000004fffffffdfffffffffffffffefffffffbffffffff0000000000000003",
));

/// Multiplicative identity.
pub(super) const ONE: FieldElement = R;

pub(super) const fn add(lhs: &FieldElement, rhs: &FieldElement) -> FieldElement {
    let a = u256_to_u64x4(lhs.0);
    let b = u256_to_u64x4(rhs.0);

    // Bit 256 of p is set, so addition can result in five words.
    let (w0, carry) = adc(a[0], b[0], 0);
    let (w1, carry) = adc(a[1], b[1], carry);
    let (w2, carry) = adc(a[2], b[2], carry);
    let (w3, w4) = adc(a[3], b[3], carry);

    // Attempt to subtract the modulus, to ensure the result is in the field.
    let modulus = u256_to_u64x4(MODULUS.0);
    let (result, _) = sub_inner(
        w0, w1, w2, w3, w4, modulus[0], modulus[1], modulus[2], modulus[3], 0,
    );
    result
}

/// Returns self * rhs mod p
pub(super) const fn multiply(lhs: &FieldElement, rhs: &FieldElement) -> FieldElement {
    // Schoolbook multiplication.
    let a = u256_to_u64x4(lhs.0);
    let b = u256_to_u64x4(rhs.0);

    let (w0, carry) = mac(0, a[0], b[0], 0);
    let (w1, carry) = mac(0, a[0], b[1], carry);
    let (w2, carry) = mac(0, a[0], b[2], carry);
    let (w3, w4) = mac(0, a[0], b[3], carry);

    let (w1, carry) = mac(w1, a[1], b[0], 0);
    let (w2, carry) = mac(w2, a[1], b[1], carry);
    let (w3, carry) = mac(w3, a[1], b[2], carry);
    let (w4, w5) = mac(w4, a[1], b[3], carry);

    let (w2, carry) = mac(w2, a[2], b[0], 0);
    let (w3, carry) = mac(w3, a[2], b[1], carry);
    let (w4, carry) = mac(w4, a[2], b[2], carry);
    let (w5, w6) = mac(w5, a[2], b[3], carry);

    let (w3, carry) = mac(w3, a[3], b[0], 0);
    let (w4, carry) = mac(w4, a[3], b[1], carry);
    let (w5, carry) = mac(w5, a[3], b[2], carry);
    let (w6, w7) = mac(w6, a[3], b[3], carry);

    montgomery_reduce(w0, w1, w2, w3, w4, w5, w6, w7)
}

// TODO(victor): Drop this function?
pub(super) fn mul_single(lhs: &FieldElement, rhs: u32) -> FieldElement {
    multiply(lhs, &FieldElement::from_u64(rhs as u64))
}

// TODO(victor): Drop this function?
pub(super) fn double(x: &FieldElement) -> FieldElement {
    add(x, x)
}

pub(super) const fn sub(lhs: &FieldElement, rhs: &FieldElement) -> FieldElement {
    let a = u256_to_u64x4(lhs.0);
    let b = u256_to_u64x4(rhs.0);
    sub_inner(a[0], a[1], a[2], a[3], 0, b[0], b[1], b[2], b[3], 0).0
}

fn from_bytes_wide(bytes: [u8; 64]) -> FieldElement {
    #[allow(clippy::unwrap_used)]
    montgomery_reduce(
        u64::from_be_bytes(bytes[0..8].try_into().unwrap()),
        u64::from_be_bytes(bytes[8..16].try_into().unwrap()),
        u64::from_be_bytes(bytes[16..24].try_into().unwrap()),
        u64::from_be_bytes(bytes[24..32].try_into().unwrap()),
        u64::from_be_bytes(bytes[32..40].try_into().unwrap()),
        u64::from_be_bytes(bytes[40..48].try_into().unwrap()),
        u64::from_be_bytes(bytes[48..56].try_into().unwrap()),
        u64::from_be_bytes(bytes[56..64].try_into().unwrap()),
    )
}

#[inline]
#[allow(clippy::too_many_arguments)]
const fn sub_inner(
    l0: u64,
    l1: u64,
    l2: u64,
    l3: u64,
    l4: u64,
    r0: u64,
    r1: u64,
    r2: u64,
    r3: u64,
    r4: u64,
) -> (FieldElement, u64) {
    let (w0, borrow) = sbb(l0, r0, 0);
    let (w1, borrow) = sbb(l1, r1, borrow);
    let (w2, borrow) = sbb(l2, r2, borrow);
    let (w3, borrow) = sbb(l3, r3, borrow);
    let (_, borrow) = sbb(l4, r4, borrow);

    // If underflow occurred on the final limb, borrow = 0xfff...fff, otherwise
    // borrow = 0x000...000. Thus, we use it as a mask to conditionally add the
    // modulus.
    let modulus = u256_to_u64x4(MODULUS.0);
    let (w0, carry) = adc(w0, modulus[0] & borrow, 0);
    let (w1, carry) = adc(w1, modulus[1] & borrow, carry);
    let (w2, carry) = adc(w2, modulus[2] & borrow, carry);
    let (w3, _) = adc(w3, modulus[3] & borrow, carry);

    (FieldElement(u64x4_to_u256([w0, w1, w2, w3])), borrow)
}

/// Montgomery Reduction
///
/// The general algorithm is:
/// ```text
/// A <- input (2n b-limbs)
/// for i in 0..n {
///     k <- A[i] p' mod b
///     A <- A + k p b^i
/// }
/// A <- A / b^n
/// if A >= p {
///     A <- A - p
/// }
/// ```
///
/// For secp256r1, we have the following simplifications:
///
/// - `p'` is 1, so our multiplicand is simply the first limb of the intermediate A.
///
/// - The first limb of p is 2^64 - 1; multiplications by this limb can be simplified
///   to a shift and subtraction:
///   ```text
///       a_i * (2^64 - 1) = a_i * 2^64 - a_i = (a_i << 64) - a_i
///   ```
///   However, because `p' = 1`, the first limb of p is multiplied by limb i of the
///   intermediate A and then immediately added to that same limb, so we simply
///   initialize the carry to limb i of the intermediate.
///
/// - The third limb of p is zero, so we can ignore any multiplications by it and just
///   add the carry.
///
/// References:
/// - Handbook of Applied Cryptography, Chapter 14
///   Algorithm 14.32
///   http://cacr.uwaterloo.ca/hac/about/chap14.pdf
///
/// - Efficient and Secure Elliptic Curve Cryptography Implementation of Curve P-256
///   Algorithm 7) Montgomery Word-by-Word Reduction
///   https://csrc.nist.gov/csrc/media/events/workshop-on-elliptic-curve-cryptography-standards/documents/papers/session6-adalier-mehmet.pdf
#[inline]
#[allow(clippy::too_many_arguments)]
const fn montgomery_reduce(
    r0: u64,
    r1: u64,
    r2: u64,
    r3: u64,
    r4: u64,
    r5: u64,
    r6: u64,
    r7: u64,
) -> FieldElement {
    let modulus = u256_to_u64x4(MODULUS.0);

    let (r1, carry) = mac(r1, r0, modulus[1], r0);
    let (r2, carry) = adc(r2, 0, carry);
    let (r3, carry) = mac(r3, r0, modulus[3], carry);
    let (r4, carry2) = adc(r4, 0, carry);

    let (r2, carry) = mac(r2, r1, modulus[1], r1);
    let (r3, carry) = adc(r3, 0, carry);
    let (r4, carry) = mac(r4, r1, modulus[3], carry);
    let (r5, carry2) = adc(r5, carry2, carry);

    let (r3, carry) = mac(r3, r2, modulus[1], r2);
    let (r4, carry) = adc(r4, 0, carry);
    let (r5, carry) = mac(r5, r2, modulus[3], carry);
    let (r6, carry2) = adc(r6, carry2, carry);

    let (r4, carry) = mac(r4, r3, modulus[1], r3);
    let (r5, carry) = adc(r5, 0, carry);
    let (r6, carry) = mac(r6, r3, modulus[3], carry);
    let (r7, r8) = adc(r7, carry2, carry);

    // Result may be within MODULUS of the correct value
    let (result, _) = sub_inner(
        r4, r5, r6, r7, r8, modulus[0], modulus[1], modulus[2], modulus[3], 0,
    );
    result
}

/// Translate a field element out of the Montgomery domain.
#[inline]
pub(super) const fn to_canonical(x: FieldElement) -> FieldElement {
    let w = u256_to_u64x4(x.0);
    montgomery_reduce(w[0], w[1], w[2], w[3], 0, 0, 0, 0)
}

/// Translate a field element into the Montgomery domain.
#[inline]
pub(super) fn to_montgomery(x: FieldElement) -> FieldElement {
    FieldElement::multiply(&x, &R2)
}

pub(super) fn to_bytes(x: FieldElement) -> FieldBytes {
    to_canonical(x).0.to_be_byte_array()
}

pub(super) fn from_uint_unchecked(w: U256) -> FieldElement {
    to_montgomery(FieldElement(w))
}

pub fn is_zero(x: &FieldElement) -> Choice {
    x.ct_eq(&FieldElement::ZERO)
}

pub(super) fn random(mut rng: impl RngCore) -> FieldElement {
    // We reduce a random 512-bit value into a 256-bit field, which results in a
    // negligible bias from the uniform distribution.
    let mut buf = [0; 64];
    rng.fill_bytes(&mut buf);
    from_bytes_wide(buf)
}
