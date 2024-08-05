//! 64-bit secp256r1 field element algorithms.

use super::{MODULUS, MODULUS_HEX};
use elliptic_curve::bigint::{Limb, U256};

pub(super) const fn add(a: U256, b: U256) -> U256 {
    let a = a.as_limbs();
    let b = b.as_limbs();

    // Bit 256 of p is set, so addition can result in five words.
    let (w0, carry) = a[0].adc(b[0], Limb::ZERO);
    let (w1, carry) = a[1].adc(b[1], carry);
    let (w2, carry) = a[2].adc(b[2], carry);
    let (w3, w4) = a[3].adc(b[3], carry);

    // Attempt to subtract the modulus, to ensure the result is in the field
    let modulus = MODULUS.0.as_limbs();

    let (result, _) = sub_inner(
        [w0, w1, w2, w3, w4],
        [modulus[0], modulus[1], modulus[2], modulus[3], Limb::ZERO],
    );
    U256::new([result[0], result[1], result[2], result[3]])
}

/// Multiplies by a single-limb integer.P
/// Multiplies the magnitude by the same value.
pub(super) fn mul_single(a: U256, rhs: u32) -> U256 {
    let mut result = U256::ZERO;
    for _i in 0..rhs {
        result = add(result, a)
    }
    result
}

/// Returns self * rhs mod p
pub(super) fn mul(a: U256, b: U256) -> U256 {
    let (lo, hi) = a.mul_wide(&b);
    let (rem, _) = U256::const_rem_wide((lo, hi), &U256::from_be_hex(MODULUS_HEX));
    rem
}

pub(super) const fn sub(a: U256, b: U256) -> U256 {
    let a = a.as_limbs();
    let b = b.as_limbs();

    let (result, _) = sub_inner(
        [a[0], a[1], a[2], a[3], Limb::ZERO],
        [b[0], b[1], b[2], b[3], Limb::ZERO],
    );
    U256::new([result[0], result[1], result[2], result[3]])
}

#[inline]
#[allow(clippy::too_many_arguments)]
const fn sub_inner(l: [Limb; 5], r: [Limb; 5]) -> ([Limb; 4], Limb) {
    let (w0, borrow) = l[0].sbb(r[0], Limb::ZERO);
    let (w1, borrow) = l[1].sbb(r[1], borrow);
    let (w2, borrow) = l[2].sbb(r[2], borrow);
    let (w3, borrow) = l[3].sbb(r[3], borrow);
    let (_, borrow) = l[4].sbb(r[4], borrow);

    // If underflow occurred on the final limb, borrow = 0xfff...fff, otherwise
    // borrow = 0x000...000. Thus, we use it as a mask to conditionally add the
    // modulus.

    let modulus = MODULUS.0.as_limbs();

    let (w0, carry) = w0.adc(modulus[0].bitand(borrow), Limb::ZERO);
    let (w1, carry) = w1.adc(modulus[1].bitand(borrow), carry);
    let (w2, carry) = w2.adc(modulus[2].bitand(borrow), carry);
    let (w3, _) = w3.adc(modulus[3].bitand(borrow), carry);

    ([w0, w1, w2, w3], borrow)
}
