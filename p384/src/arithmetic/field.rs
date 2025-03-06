//! Field arithmetic modulo p = 2^{384} − 2^{128} − 2^{96} + 2^{32} − 1
//!
//! Arithmetic implementations have been synthesized using fiat-crypto.
//!
//! # License
//!
//! Copyright (c) 2015-2020 the fiat-crypto authors
//!
//! fiat-crypto is distributed under the terms of the MIT License, the
//! Apache License (Version 2.0), and the BSD 1-Clause License;
//! users may pick which license to apply.

#![allow(
    clippy::should_implement_trait,
    clippy::suspicious_op_assign_impl,
    clippy::unused_unit,
    clippy::unnecessary_cast,
    clippy::too_many_arguments,
    clippy::identity_op
)]

#[cfg_attr(target_pointer_width = "32", path = "field/p384_32.rs")]
#[cfg_attr(target_pointer_width = "64", path = "field/p384_64.rs")]
mod field_impl;

use self::field_impl::*;
use crate::{FieldBytes, NistP384};
use core::{
    fmt::{self, Debug},
    iter::{Product, Sum},
    ops::{AddAssign, MulAssign, Neg, SubAssign},
};
use elliptic_curve::ops::Invert;
use elliptic_curve::{
    bigint::{Limb, U384},
    ff::PrimeField,
    subtle::{Choice, ConstantTimeEq, CtOption},
};
use primeorder::impl_bernstein_yang_invert;

/// Constant representing the modulus
/// p = 2^{384} − 2^{128} − 2^{96} + 2^{32} − 1
pub(crate) const MODULUS: U384 = U384::from_be_hex(FieldElement::MODULUS);

/// Element of the secp384r1 base field used for curve coordinates.
#[derive(Clone, Copy)]
pub struct FieldElement(pub(super) U384);

primeorder::impl_mont_field_element!(
    NistP384,
    FieldElement,
    FieldBytes,
    U384,
    MODULUS,
    fiat_p384_montgomery_domain_field_element,
    fiat_p384_from_montgomery,
    fiat_p384_to_montgomery,
    fiat_p384_add,
    fiat_p384_sub,
    fiat_p384_mul,
    fiat_p384_opp,
    fiat_p384_square
);

#[cfg(all(target_os = "zkvm", target_arch = "riscv32"))]
use primeorder::__risc0::FieldElement384;

#[cfg(all(target_os = "zkvm", target_arch = "riscv32"))]
use risc0_bigint2::field::FIELD_384_WIDTH_WORDS;

#[cfg(all(target_os = "zkvm", target_arch = "riscv32"))]
const R_2_LE: FieldElement384<NistP384> = FieldElement384::new_unchecked([
    0x00000001, 0xfffffffe, 0x00000000, 0x00000002, 0x00000000, 0xfffffffe, 0x00000000, 0x00000002,
    0x00000001, 0x00000000, 0x00000000, 0x00000000,
]);

impl FieldElement {
    #[cfg(all(target_os = "zkvm", target_arch = "riscv32"))]
    pub(crate) fn from_words_le(fe: [u32; FIELD_384_WIDTH_WORDS]) -> Self {
        let fe = FieldElement384::new_unchecked(fe);

        // Convert to montgomery form with aR mod p
        let mut mont = FieldElement384::default();

        // This mul will check if the result is within the modulus.
        fe.mul(&R_2_LE, &mut mont);

        let uint = U384::from_le_slice(bytemuck::cast_slice::<u32, u8>(&mont.data));

        Self(uint)
    }

    #[cfg(all(target_os = "zkvm", target_arch = "riscv32"))]
    pub(crate) fn to_words_le(&self) -> [u32; FIELD_384_WIDTH_WORDS] {
        use crate::elliptic_curve::bigint::Encoding;
        // NOTE: this from mont conversion could be accelerated, but it's very little cycles.
        let canonical = self.to_canonical();
        let input = canonical.to_le_bytes();
        let array = bytemuck::cast::<_, [u32; 12]>(input);

        array
    }

    /// Compute [`FieldElement`] inversion: `1 / self`.
    pub fn invert(&self) -> CtOption<Self> {
        #[cfg(all(target_os = "zkvm", target_arch = "riscv32"))]
        {
            use crate::elliptic_curve::bigint::Encoding;

            // NOTE: This is not a constant time operation, as inverting zero in the zkvm is not
            // possible as it will panic in the host.
            if self.is_zero().into() {
                return CtOption::new(FieldElement::ZERO, Choice::from(0));
            } else {
                let input_words = self.to_words_le();
                let mut output = [0u32; FIELD_384_WIDTH_WORDS];
                risc0_bigint2::field::modinv_384(
                    &input_words,
                    &crate::__risc0::SECP384R1_PRIME,
                    &mut output,
                );
                let element = FieldElement::from_words_le(output);
                return CtOption::new(element, Choice::from(1));
            }
        }
        
        #[cfg(not(all(target_os = "zkvm", target_arch = "riscv32")))]
        CtOption::new(self.invert_unchecked(), !self.is_zero())
    }

    /// Returns the multiplicative inverse of self.
    ///
    /// Does not check that self is non-zero.
    const fn invert_unchecked(&self) -> Self {
        let words = impl_bernstein_yang_invert!(
            self.0.as_words(),
            Self::ONE.0.to_words(),
            384,
            U384::LIMBS,
            Limb,
            fiat_p384_from_montgomery,
            fiat_p384_mul,
            fiat_p384_opp,
            fiat_p384_divstep_precomp,
            fiat_p384_divstep,
            fiat_p384_msat,
            fiat_p384_selectznz,
        );

        Self(U384::from_words(words))
    }

    /// Returns the square root of self mod p, or `None` if no square root
    /// exists.
    pub fn sqrt(&self) -> CtOption<Self> {
        // p mod 4 = 3 -> compute sqrt(x) using x^((p+1)/4) =
        // x^9850501549098619803069760025035903451269934817616361666987073351061430442874217582261816522064734500465401743278080
        let t1 = *self;
        let t10 = t1.square();
        let t11 = t1 * t10;
        let t110 = t11.square();
        let t111 = t1 * t110;
        let t111000 = t111.sqn(3);
        let t111111 = t111 * t111000;
        let t1111110 = t111111.square();
        let t1111111 = t1 * t1111110;
        let x12 = t1111110.sqn(5) * t111111;
        let x24 = x12.sqn(12) * x12;
        let x31 = x24.sqn(7) * t1111111;
        let x32 = x31.square() * t1;
        let x63 = x32.sqn(31) * x31;
        let x126 = x63.sqn(63) * x63;
        let x252 = x126.sqn(126) * x126;
        let x255 = x252.sqn(3) * t111;
        let x = ((x255.sqn(33) * x32).sqn(64) * t1).sqn(30);
        CtOption::new(x, x.square().ct_eq(&t1))
    }

    /// Returns self^(2^n) mod p.
    fn sqn(&self, n: usize) -> Self {
        let mut x = *self;
        for _ in 0..n {
            x = x.square();
        }
        x
    }
}

impl PrimeField for FieldElement {
    type Repr = FieldBytes;

    const MODULUS: &'static str = "fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffeffffffff0000000000000000ffffffff";
    const NUM_BITS: u32 = 384;
    const CAPACITY: u32 = 383;
    const TWO_INV: Self = Self::from_u64(2).invert_unchecked();
    const MULTIPLICATIVE_GENERATOR: Self = Self::from_u64(19);
    const S: u32 = 1;
    const ROOT_OF_UNITY: Self = Self::from_hex("fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffeffffffff0000000000000000fffffffe");
    const ROOT_OF_UNITY_INV: Self = Self::ROOT_OF_UNITY.invert_unchecked();
    const DELTA: Self = Self::from_u64(49);

    #[inline]
    fn from_repr(bytes: FieldBytes) -> CtOption<Self> {
        Self::from_bytes(&bytes)
    }

    #[inline]
    fn to_repr(&self) -> FieldBytes {
        self.to_bytes()
    }

    #[inline]
    fn is_odd(&self) -> Choice {
        self.is_odd()
    }
}

impl Debug for FieldElement {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "FieldElement(0x{:X})", &self.0)
    }
}

impl Invert for FieldElement {
    type Output = CtOption<Self>;

    fn invert(&self) -> CtOption<Self> {
        self.invert()
    }
}

#[cfg(test)]
mod tests {
    use super::FieldElement;
    use elliptic_curve::ff::PrimeField;
    use primeorder::{
        impl_field_identity_tests, impl_field_invert_tests, impl_field_sqrt_tests,
        impl_primefield_tests,
    };

    /// t = (modulus - 1) >> S
    const T: [u64; 6] = [
        0x000000007fffffff,
        0x7fffffff80000000,
        0xffffffffffffffff,
        0xffffffffffffffff,
        0xffffffffffffffff,
        0x7fffffffffffffff,
    ];

    impl_field_identity_tests!(FieldElement);
    impl_field_invert_tests!(FieldElement);
    impl_field_sqrt_tests!(FieldElement);
    impl_primefield_tests!(FieldElement, T);
}
