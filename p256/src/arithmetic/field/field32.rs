//! 32-bit secp256r1 field element algorithms.

use super::{MODULUS, MODULUS_HEX};
use elliptic_curve::bigint::{Limb, U256};

pub(super) const fn add(a: U256, b: U256) -> U256 {
    let a = a.as_limbs();
    let b = b.as_limbs();

    // Bit 256 of p is set, so addition can result in nine words.
    // let (w0, carry) = adc(a[0], b[0], 0);
    let (w0, carry) = a[0].adc(b[0], Limb::ZERO);
    let (w1, carry) = a[1].adc(b[1], carry);
    let (w2, carry) = a[2].adc(b[2], carry);
    let (w3, carry) = a[3].adc(b[3], carry);
    let (w4, carry) = a[4].adc(b[4], carry);
    let (w5, carry) = a[5].adc(b[5], carry);
    let (w6, carry) = a[6].adc(b[6], carry);
    let (w7, w8) = a[7].adc(b[7], carry);
    // Attempt to subtract the modulus, to ensure the result is in the field.
    let modulus = MODULUS.0.as_limbs();

    let (result, _) = sub_inner(
        [w0, w1, w2, w3, w4, w5, w6, w7, w8],
        [
            modulus[0],
            modulus[1],
            modulus[2],
            modulus[3],
            modulus[4],
            modulus[5],
            modulus[6],
            modulus[7],
            Limb::ZERO,
        ],
    );
    U256::new([
        result[0], result[1], result[2], result[3], result[4], result[5], result[6], result[7],
    ])
}

/// Multiplies by a single-limb integer.
/// Multiplies the magnitude by the same value.
pub(super) fn mul_single(a: U256, rhs: u32) -> U256 {
    let a_limbs = a.as_limbs();
    let rhs_limb = Limb::from_u32(rhs);
    let (w0, carry) = a[0].mac(b[0], Limb::ZERO);
    let (w1, carry) = a[1].mac(b[1], carry);
    let (w2, carry) = a[2].mac(b[2], carry);
    let (w3, carry) = a[3].mac(b[3], carry);
    let (w4, carry) = a[4].mac(b[4], carry);
    let (w5, carry) = a[5].mac(b[5], carry);
    let (w6, carry) = a[6].mac(b[6], carry);
    let (w7, w8) = a[7].mac(b[7], carry);

    // Reduce the carry mod prime
    let carry = U256::from(w8);
    let (reduced_carry, _) = carry.const_rem(&MODULUS.0);

    // Modular addition of non-carry and reduced carry
    let non_carries = U256::new([w0, w1, w2, w3, w4, w5, w6, w7]);
    add(non_carries, reduced_carry)
}

/// Returns self * rhs mod p
pub(super) const fn mul(a: U256, b: U256) -> U256 {
    let (lo, hi): (U256, U256) = a.mul_wide(&b);
    let (rem, _) = U256::const_rem_wide((lo, hi), &U256::from_be_hex(MODULUS_HEX));
    rem
}

pub(super) const fn sub(a: U256, b: U256) -> U256 {
    let a = a.as_limbs();
    let b = b.as_limbs();

    let (result, _) = sub_inner(
        [a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], Limb::ZERO],
        [b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7], Limb::ZERO],
    );
    U256::new([
        result[0], result[1], result[2], result[3], result[4], result[5], result[6], result[7],
    ])
}

#[inline]
#[allow(clippy::too_many_arguments)]
const fn sub_inner(l: [Limb; 9], r: [Limb; 9]) -> ([Limb; 8], Limb) {
    let (w0, borrow) = l[0].sbb(r[0], Limb::ZERO);
    let (w1, borrow) = l[1].sbb(r[1], borrow);
    let (w2, borrow) = l[2].sbb(r[2], borrow);
    let (w3, borrow) = l[3].sbb(r[3], borrow);
    let (w4, borrow) = l[4].sbb(r[4], borrow);
    let (w5, borrow) = l[5].sbb(r[5], borrow);
    let (w6, borrow) = l[6].sbb(r[6], borrow);
    let (w7, borrow) = l[7].sbb(r[7], borrow);
    let (_, borrow) = l[8].sbb(r[8], borrow);

    // If underflow occurred on the final limb, borrow = 0xfff...fff, otherwise
    // borrow = 0x000...000. Thus, we use it as a mask to conditionally add
    // the modulus.

    let modulus = MODULUS.0.as_limbs();

    let (w0, carry) = w0.adc(modulus[0].bitand(borrow), Limb::ZERO);
    let (w1, carry) = w1.adc(modulus[1].bitand(borrow), carry);
    let (w2, carry) = w2.adc(modulus[2].bitand(borrow), carry);
    let (w3, carry) = w3.adc(modulus[3].bitand(borrow), carry);
    let (w4, carry) = w4.adc(modulus[4].bitand(borrow), carry);
    let (w5, carry) = w5.adc(modulus[5].bitand(borrow), carry);
    let (w6, carry) = w6.adc(modulus[6].bitand(borrow), carry);
    let (w7, _) = w7.adc(modulus[7].bitand(borrow), carry);

    ([w0, w1, w2, w3, w4, w5, w6, w7], borrow)
}
