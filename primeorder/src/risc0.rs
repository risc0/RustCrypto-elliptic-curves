use crate::FieldBytes;
use crate::{affine::AffinePoint, projective::ProjectivePoint};
use core::marker::PhantomData;
use core::ops::Deref;
use elliptic_curve::consts::U32;
use elliptic_curve::generic_array::GenericArray;
use elliptic_curve::subtle::{Choice, ConditionallySelectable};
use elliptic_curve::{PrimeField, Scalar};
use risc0_bigint2::ec;

use crate::PrimeCurveParams;
use crate::PrimeCurveParams256;
use crate::PrimeCurveParams384;

/// Representation of a field element in raw bytes form. This is not in montgomery form.
#[derive(Copy, Clone, Default, Debug, PartialEq, Eq)]
pub struct FieldElement256<C> {
    pub data: [u32; 8],
    _phantom: PhantomData<C>,
}

impl<C> Deref for FieldElement256<C> {
    type Target = [u32; 8];

    fn deref(&self) -> &[u32; 8] {
        &self.data
    }
}

impl<C: elliptic_curve::Curve> From<&FieldBytes<C>> for FieldElement256<C> {
    fn from(data: &FieldBytes<C>) -> Self {
        let mut words = [0u32; 8];

        // Process 4 bytes at a time to create little-endian u32 words
        for (i, chunk) in data.chunks(4).enumerate() {
            // Convert each big-endian chunk to a little-endian u32
            words[7 - i] = u32::from_be_bytes(chunk.try_into().unwrap());
        }

        Self::new_unchecked(words)
    }
}

impl<C: PrimeCurveParams> From<FieldElement256<C>> for GenericArray<u8, C::FieldBytesSize> {
    fn from(data: FieldElement256<C>) -> Self {
        let bytes_slice = bytemuck::cast_slice::<u32, u8>(&data.data);
        GenericArray::from_iter(bytes_slice.iter().copied().rev())
    }
}

impl<C: Copy> ConditionallySelectable for FieldElement256<C> {
    #[inline]
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        let mut output = *a;
        output.conditional_assign(b, choice);
        output
    }

    fn conditional_assign(&mut self, other: &Self, choice: Choice) {
        for (a_i, b_i) in self.data.iter_mut().zip(other.data.iter()) {
            a_i.conditional_assign(b_i, choice)
        }
    }
}

impl<C> FieldElement256<C> {
    pub const fn new_unchecked(data: [u32; 8]) -> Self {
        Self {
            data,
            _phantom: PhantomData,
        }
    }
}

impl<C> FieldElement256<C>
where
    C: PrimeCurveParams256,
{
    #[inline]
    pub fn mul_unchecked(&self, rhs: &Self, result: &mut Self) {
        risc0_bigint2::field::unchecked::modmul_256(
            &self.data,
            &rhs.data,
            &C::PRIME_LE_WORDS,
            &mut result.data,
        );
    }

    #[inline]
    pub fn mul(&self, rhs: &Self, result: &mut Self) {
        risc0_bigint2::field::modmul_256(
            &self.data,
            &rhs.data,
            &C::PRIME_LE_WORDS,
            &mut result.data,
        );
    }

    #[inline]
    pub fn add_unchecked(&self, rhs: &Self, result: &mut Self) {
        risc0_bigint2::field::unchecked::modadd_256(
            &self.data,
            &rhs.data,
            &C::PRIME_LE_WORDS,
            &mut result.data,
        );
    }

    /// Calculate the square root of the field element, using the provided buffer as scratch, and
    /// writing the result to the `result` parameter.
    pub fn sqrt_unchecked(&self, scratch: &mut Self, result: &mut Self) {
        // New buffers to keep temporary values.
        let mut scratch_1 = Self::default();
        let mut scratch_2 = Self::default();

        // let t11 = self.mul(&self.square());
        self.square(scratch);
        self.mul_unchecked(scratch, result);

        // result = t11
        // let t1111 = t11.mul(&t11.sqn(2));
        result.sqn(2, (&mut scratch_1, &mut scratch_2), scratch);
        result.mul_unchecked(scratch, &mut scratch_1);

        // scratch_1 = t1111
        // let t11111111 = t1111.mul(&t1111.sqn(4));
        scratch_1.sqn(4, (&mut scratch_2, result), scratch);
        scratch_1.mul_unchecked(scratch, result);

        // result = t11111111
        // let x16 = t11111111.sqn(8).mul(&t11111111);
        result.sqn(8, (&mut scratch_1, &mut scratch_2), scratch);
        result.mul_unchecked(scratch, &mut scratch_1);

        // scratch_1 = x16
        // let sqrt = x16
        //     .sqn(16)
        scratch_1.sqn(16, (&mut scratch_2, result), scratch);
        //     .mul(&x16)
        scratch.mul_unchecked(&scratch_1, result);
        //     .sqn(32)
        result.sqn(32, (&mut scratch_1, &mut scratch_2), scratch);
        //     .mul(self)
        scratch.mul_unchecked(self, &mut scratch_1);
        //     .sqn(96)
        scratch_1.sqn(96, (&mut scratch_2, result), scratch);
        //     .mul(self)
        scratch.mul_unchecked(self, &mut scratch_2);
        //     .sqn(94);
        // Last result is written to the result buffer.
        scratch_2.sqn(94, (&mut scratch_1, scratch), result);
    }

    /// Returns self^(2^n) mod p.
    ///
    /// This implementation is designed to avoid any memcpy of the buffers for intermediate ops.
    fn sqn(&self, n: usize, scratch: (&mut Self, &mut Self), result: &mut Self) {
        let mut x = scratch.0;
        let mut buffer = scratch.1;

        if n == 1 {
            // self^(2^1) = self^2
            self.square(result);
            return;
        } else if n == 0 {
            // self^(1) = self
            *result = *self;
            return;
        }

        // write value to a scratch buffer.
        self.square(x);

        // Square n - 2 times.
        let mut i = 2;
        while i < n {
            x.square(buffer);
            i += 1;
            // Swap scratch buffers, to set x to the squared value.
            core::mem::swap(&mut x, &mut buffer);
        }

        // Write final square to result buffer.
        x.square(result);
    }

    /// Returns self^2 mod p
    pub fn square(&self, result: &mut Self) {
        self.mul_unchecked(self, result);
    }
}

fn bytes_to_u32_words_le(bytes: &[u8]) -> [u32; 8] {
    let mut words = [0u32; 8];

    // Process 4 bytes at a time to create little-endian u32 words
    for (i, chunk) in bytes.chunks(4).enumerate() {
        // Convert each big-endian chunk to a little-endian u32
        words[7 - i] = u32::from_be_bytes(chunk.try_into().unwrap());
    }

    words
}

pub fn felt_to_u32_words_le<C>(data: &C::FieldElement) -> [u32; 8]
where
    C: PrimeCurveParams256,
{
    bytes_to_u32_words_le(data.to_repr().as_slice())
}

#[inline]
fn affine_to_r0_affine<C>(affine: &AffinePoint<C>) -> ec::AffinePoint<8, C>
where
    C: PrimeCurveParams256,
{
    if bool::from(affine.is_identity()) {
        return ec::AffinePoint::IDENTITY;
    }

    let x = felt_to_u32_words_le::<C>(&affine.x);
    let y = felt_to_u32_words_le::<C>(&affine.y);
    ec::AffinePoint::new_unchecked(x, y)
}

pub fn projective_to_affine<C>(p: &ProjectivePoint<C>) -> ec::AffinePoint<8, C>
where
    C: PrimeCurveParams256,
{
    let aff = p.to_affine();
    affine_to_r0_affine(&aff)
}

/// Public helper for curves implementing zkvm_to_affine (256-bit).
/// Converts projective point to affine using accelerated field inversion.
pub fn zkvm_to_affine_impl_256<C>(point: &ProjectivePoint<C>) -> AffinePoint<C>
where
    C: PrimeCurveParams256,
{
    use elliptic_curve::{Field, PrimeField};

    if point.z.is_zero().into() {
        return AffinePoint::IDENTITY;
    }
    let z = felt_to_u32_words_le::<C>(&point.z);
    let mut z_inv = [0u32; 8];
    risc0_bigint2::field::unchecked::modinv_256(&z, &C::PRIME_LE_WORDS, &mut z_inv);

    let mut buffer = [0u32; 8];
    let x_buffer = felt_to_u32_words_le::<C>(&point.x);
    let y_buffer = felt_to_u32_words_le::<C>(&point.y);
    risc0_bigint2::field::unchecked::modmul_256(
        &x_buffer,
        &z_inv,
        &C::PRIME_LE_WORDS,
        &mut buffer,
    );

    let x = C::from_u32_words_le(buffer);

    risc0_bigint2::field::unchecked::modmul_256(
        &y_buffer,
        &z_inv,
        &C::PRIME_LE_WORDS,
        &mut buffer,
    );
    let y = C::from_u32_words_le(buffer);
    AffinePoint { x, y, infinity: 0 }
}

/// Public helper for curves implementing zkvm_decompress (256-bit).
/// Decompresses a point from its x-coordinate using accelerated field operations.
pub fn zkvm_decompress_impl_256<C>(
    x_bytes: &FieldBytes<C>,
    y_is_odd: elliptic_curve::subtle::Choice,
) -> elliptic_curve::subtle::CtOption<AffinePoint<C>>
where
    C: PrimeCurveParams256 + elliptic_curve::Curve<FieldBytesSize = elliptic_curve::consts::U32>,
{
    use elliptic_curve::subtle::{ConditionallySelectable, ConstantTimeEq, CtOption};

    // Note: buffers are kept separate for each OP as the result pointer cannot equal one
    // of the input pointers.
    let mut scratch = FieldElement256::<C>::default();
    let mut acc = FieldElement256::<C>::default();
    let mut scratch_1 = FieldElement256::<C>::from(x_bytes);

    // x checked to be in the field.
    C::FieldElement::from_repr(*x_bytes).and_then(|x| {
        // x * &x * &x
        scratch_1.mul_unchecked(&scratch_1, &mut scratch);
        scratch.mul_unchecked(&scratch_1, &mut acc);

        // + &(C::EQUATION_A * &x)
        scratch_1.mul_unchecked(&C::EQUATION_A_LE, &mut scratch);
        // Can re-use x as a buffer, no longer needed.
        scratch.add_unchecked(&acc, &mut scratch_1);

        // + &C::EQUATION_B
        scratch_1.add_unchecked(&C::EQUATION_B_LE, &mut scratch);

        // Sqrt implementation. Not separated into another function to allow
        // re-using buffers.
        scratch.sqrt_unchecked(&mut scratch_1, &mut acc);

        // Check that the square root is correct.
        acc.square(&mut scratch_1);

        let sqrt = CtOption::new(acc, elliptic_curve::subtle::Choice::from(scratch_1.eq(&scratch) as u8));

        // Checked that the result is within the field.
        sqrt.and_then(|sqrt| {
            C::FieldElement::from_repr(sqrt.into()).map(|beta| {
                let y = C::FieldElement::conditional_select(
                    &-beta,
                    &beta,
                    beta.is_odd().ct_eq(&y_is_odd),
                );

                AffinePoint { x, y, infinity: 0 }
            })
        })
    })
}

/// Public helper for curves implementing zkvm_to_affine (384-bit).
/// Converts projective point to affine using accelerated field inversion.
pub fn zkvm_to_affine_impl_384<C>(point: &ProjectivePoint<C>) -> AffinePoint<C>
where
    C: PrimeCurveParams384,
{
    use elliptic_curve::{Field, PrimeField};

    if point.z.is_zero().into() {
        return AffinePoint::IDENTITY;
    }
    let z = felt_to_u32_words_le_384::<C>(&point.z);
    let mut z_inv = [0u32; 12];
    risc0_bigint2::field::unchecked::modinv_384(&z, &C::PRIME_LE_WORDS, &mut z_inv);

    let mut buffer = [0u32; 12];
    let x_buffer = felt_to_u32_words_le_384::<C>(&point.x);
    let y_buffer = felt_to_u32_words_le_384::<C>(&point.y);
    risc0_bigint2::field::unchecked::modmul_384(
        &x_buffer,
        &z_inv,
        &C::PRIME_LE_WORDS,
        &mut buffer,
    );

    let x = C::from_u32_words_le(buffer);

    risc0_bigint2::field::unchecked::modmul_384(
        &y_buffer,
        &z_inv,
        &C::PRIME_LE_WORDS,
        &mut buffer,
    );
    let y = C::from_u32_words_le(buffer);
    AffinePoint { x, y, infinity: 0 }
}

/// Public helper for curves implementing zkvm_decompress (384-bit).
/// Decompresses a point from its x-coordinate using accelerated field operations.
pub fn zkvm_decompress_impl_384<C>(
    x_bytes: &FieldBytes<C>,
    y_is_odd: elliptic_curve::subtle::Choice,
) -> elliptic_curve::subtle::CtOption<AffinePoint<C>>
where
    C: PrimeCurveParams384 + elliptic_curve::Curve<FieldBytesSize = elliptic_curve::consts::U48>,
{
    use elliptic_curve::subtle::{ConditionallySelectable, ConstantTimeEq, CtOption};

    // Note: buffers are kept separate for each OP as the result pointer cannot equal one
    // of the input pointers.
    let mut scratch = FieldElement384::<C>::default();
    let mut acc = FieldElement384::<C>::default();
    let mut scratch_1 = FieldElement384::<C>::from(x_bytes);

    // x checked to be in the field.
    C::FieldElement::from_repr(*x_bytes).and_then(|x| {
        // x * &x * &x
        scratch_1.mul_unchecked(&scratch_1, &mut scratch);
        scratch.mul_unchecked(&scratch_1, &mut acc);

        // + &(C::EQUATION_A * &x)
        scratch_1.mul_unchecked(&C::EQUATION_A_LE, &mut scratch);
        // Can re-use x as a buffer, no longer needed.
        scratch.add_unchecked(&acc, &mut scratch_1);

        // + &C::EQUATION_B
        scratch_1.add_unchecked(&C::EQUATION_B_LE, &mut scratch);

        // Sqrt implementation.
        scratch.sqrt_unchecked(&mut scratch_1, &mut acc);

        // Check that the square root is correct.
        acc.square(&mut scratch_1);

        let sqrt = CtOption::new(acc, elliptic_curve::subtle::Choice::from(scratch_1.eq(&scratch) as u8));

        // Checked that the result is within the field.
        sqrt.and_then(|sqrt| {
            C::FieldElement::from_repr(sqrt.into()).map(|beta| {
                let y = C::FieldElement::conditional_select(
                    &-beta,
                    &beta,
                    beta.is_odd().ct_eq(&y_is_odd),
                );

                AffinePoint { x, y, infinity: 0 }
            })
        })
    })
}

pub fn affine_to_projective<C>(affine: &ec::AffinePoint<8, C>) -> ProjectivePoint<C>
where
    C: PrimeCurveParams256,
{
    if let Some(value) = affine.as_u32s() {
        // This should only not be within the modulus with a malicious host, panic in that case.
        let x = C::from_u32_words_le(value[0]);
        let y = C::from_u32_words_le(value[1]);

        let affine = AffinePoint { x, y, infinity: 0 };
        ProjectivePoint::from(affine)
    } else {
        ProjectivePoint::IDENTITY
    }
}

pub fn scalar_to_words<C>(s: &Scalar<C>) -> [u32; 8]
where
    C: PrimeCurveParams256,
{
    bytes_to_u32_words_le(s.to_repr().as_slice())
}

pub mod ec_impl {
    use super::*;

    pub fn mul<C>(lhs: &ProjectivePoint<C>, rhs: &Scalar<C>) -> ProjectivePoint<C>
    where
        C: PrimeCurveParams256,
    {
        let scalar = scalar_to_words::<C>(rhs);
        let affine = projective_to_affine::<C>(lhs);

        let mut result = ec::AffinePoint::new_unchecked([0u32; 8], [0u32; 8]);
        affine.mul(&scalar, &mut result);
        return affine_to_projective(&result);
    }

    pub fn add<C>(lhs: &ProjectivePoint<C>, rhs: &ProjectivePoint<C>) -> ProjectivePoint<C>
    where
        C: PrimeCurveParams256,
    {
        let lhs = projective_to_affine::<C>(lhs);
        let rhs = projective_to_affine::<C>(rhs);

        let mut result = ec::AffinePoint::new_unchecked([0u32; 8], [0u32; 8]);
        lhs.add(&rhs, &mut result);
        return affine_to_projective(&result);
    }

    #[inline]
    pub fn add_mixed<C>(lhs: &ProjectivePoint<C>, rhs: &AffinePoint<C>) -> ProjectivePoint<C>
    where
        C: PrimeCurveParams256,
    {
        let lhs = projective_to_affine::<C>(lhs);
        let rhs = affine_to_r0_affine(rhs);

        let mut result = ec::AffinePoint::new_unchecked([0u32; 8], [0u32; 8]);
        lhs.add(&rhs, &mut result);
        return affine_to_projective(&result);
    }

    pub fn double<C>(point: &ProjectivePoint<C>) -> ProjectivePoint<C>
    where
        C: PrimeCurveParams256,
    {
        let point = projective_to_affine::<C>(point);

        let mut result = ec::AffinePoint::new_unchecked([0u32; 8], [0u32; 8]);
        point.double(&mut result);
        return affine_to_projective(&result);
    }
}

// ============================================================================
// 384-bit curve support
// ============================================================================

use elliptic_curve::consts::U48;

/// Representation of a 384-bit field element in raw bytes form. This is not in montgomery form.
#[derive(Copy, Clone, Default, Debug, PartialEq, Eq)]
pub struct FieldElement384<C> {
    pub data: [u32; 12],
    _phantom: PhantomData<C>,
}

impl<C> Deref for FieldElement384<C> {
    type Target = [u32; 12];

    fn deref(&self) -> &[u32; 12] {
        &self.data
    }
}

impl<C: elliptic_curve::Curve<FieldBytesSize = U48>> From<&FieldBytes<C>> for FieldElement384<C> {
    fn from(data: &FieldBytes<C>) -> Self {
        let mut words = [0u32; 12];

        // Process 4 bytes at a time to create little-endian u32 words
        for (i, chunk) in data.chunks(4).enumerate() {
            // Convert each big-endian chunk to a little-endian u32
            words[11 - i] = u32::from_be_bytes(chunk.try_into().unwrap());
        }

        Self::new_unchecked(words)
    }
}

impl<C: PrimeCurveParams384> From<FieldElement384<C>> for GenericArray<u8, C::FieldBytesSize> {
    fn from(data: FieldElement384<C>) -> Self {
        let bytes_slice = bytemuck::cast_slice::<u32, u8>(&data.data);
        GenericArray::from_iter(bytes_slice.iter().copied().rev())
    }
}

impl<C: Copy> ConditionallySelectable for FieldElement384<C> {
    #[inline]
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        let mut output = *a;
        output.conditional_assign(b, choice);
        output
    }

    fn conditional_assign(&mut self, other: &Self, choice: Choice) {
        for (a_i, b_i) in self.data.iter_mut().zip(other.data.iter()) {
            a_i.conditional_assign(b_i, choice)
        }
    }
}

impl<C> FieldElement384<C> {
    pub const fn new_unchecked(data: [u32; 12]) -> Self {
        Self {
            data,
            _phantom: PhantomData,
        }
    }
}

impl<C> FieldElement384<C>
where
    C: PrimeCurveParams384,
{
    #[inline]
    pub fn mul_unchecked(&self, rhs: &Self, result: &mut Self) {
        risc0_bigint2::field::unchecked::modmul_384(
            &self.data,
            &rhs.data,
            &C::PRIME_LE_WORDS,
            &mut result.data,
        );
    }

    #[inline]
    pub fn mul(&self, rhs: &Self, result: &mut Self) {
        risc0_bigint2::field::modmul_384(
            &self.data,
            &rhs.data,
            &C::PRIME_LE_WORDS,
            &mut result.data,
        );
    }

    #[inline]
    pub fn add_unchecked(&self, rhs: &Self, result: &mut Self) {
        risc0_bigint2::field::unchecked::modadd_384(
            &self.data,
            &rhs.data,
            &C::PRIME_LE_WORDS,
            &mut result.data,
        );
    }

    /// Calculate the square root of the field element for P-384.
    /// Uses x^((p+1)/4) since p ≡ 3 (mod 4).
    /// Based on the addition chain from p384/src/arithmetic/field.rs
    pub fn sqrt_unchecked(&self, scratch: &mut Self, result: &mut Self) {
        // We need multiple buffers for intermediate values
        let mut buf_a = Self::default();
        let mut buf_b = Self::default();
        let mut buf_c = Self::default();
        let mut buf_d = Self::default();

        // t1 = self (input)
        // t10 = t1^2
        self.square(&mut buf_a);  // buf_a = t10

        // t11 = t1 * t10 = self^3
        self.mul_unchecked(&buf_a, &mut buf_b);  // buf_b = t11

        // t110 = t11^2 = self^6
        buf_b.square(&mut buf_c);  // buf_c = t110

        // t111 = t1 * t110 = self^7
        self.mul_unchecked(&buf_c, &mut buf_a);  // buf_a = t111
        let t111 = buf_a;  // save t111

        // t111000 = t111^(2^3) = self^56
        buf_a.sqn(3, (&mut buf_b, &mut buf_c), &mut buf_d);  // buf_d = t111000

        // t111111 = t111 * t111000 = self^63
        t111.mul_unchecked(&buf_d, &mut buf_a);  // buf_a = t111111
        let t111111 = buf_a;  // save t111111

        // t1111110 = t111111^2 = self^126
        buf_a.square(&mut buf_b);  // buf_b = t1111110

        // t1111111 = t1 * t1111110 = self^127
        self.mul_unchecked(&buf_b, &mut buf_c);  // buf_c = t1111111
        let t1111111 = buf_c;  // save t1111111

        // x12 = t1111110^(2^5) * t111111
        buf_b.sqn(5, (&mut buf_a, &mut buf_d), &mut buf_c);  // buf_c = t1111110^32
        buf_c.mul_unchecked(&t111111, &mut buf_a);  // buf_a = x12
        let x12 = buf_a;

        // x24 = x12^(2^12) * x12
        buf_a.sqn(12, (&mut buf_b, &mut buf_c), &mut buf_d);  // buf_d = x12^(2^12)
        buf_d.mul_unchecked(&x12, &mut buf_a);  // buf_a = x24

        // x31 = x24^(2^7) * t1111111
        buf_a.sqn(7, (&mut buf_b, &mut buf_c), &mut buf_d);  // buf_d = x24^(2^7)
        buf_d.mul_unchecked(&t1111111, &mut buf_a);  // buf_a = x31
        let x31 = buf_a;

        // x32 = x31^2 * t1
        buf_a.square(&mut buf_b);  // buf_b = x31^2
        self.mul_unchecked(&buf_b, &mut buf_c);  // buf_c = x32
        let x32 = buf_c;

        // x63 = x32^(2^31) * x31
        buf_c.sqn(31, (&mut buf_a, &mut buf_b), &mut buf_d);  // buf_d = x32^(2^31)
        buf_d.mul_unchecked(&x31, &mut buf_a);  // buf_a = x63
        let x63 = buf_a;

        // x126 = x63^(2^63) * x63
        buf_a.sqn(63, (&mut buf_b, &mut buf_c), &mut buf_d);  // buf_d = x63^(2^63)
        buf_d.mul_unchecked(&x63, &mut buf_a);  // buf_a = x126
        let x126 = buf_a;

        // x252 = x126^(2^126) * x126
        buf_a.sqn(126, (&mut buf_b, &mut buf_c), &mut buf_d);  // buf_d = x126^(2^126)
        buf_d.mul_unchecked(&x126, &mut buf_a);  // buf_a = x252

        // x255 = x252^(2^3) * t111
        buf_a.sqn(3, (&mut buf_b, &mut buf_c), &mut buf_d);  // buf_d = x252^(2^3)
        buf_d.mul_unchecked(&t111, &mut buf_a);  // buf_a = x255

        // x = ((x255^(2^33) * x32)^(2^64) * t1)^(2^30)
        buf_a.sqn(33, (&mut buf_b, &mut buf_c), &mut buf_d);  // buf_d = x255^(2^33)
        buf_d.mul_unchecked(&x32, &mut buf_a);  // buf_a = x255^(2^33) * x32
        buf_a.sqn(64, (&mut buf_b, &mut buf_c), &mut buf_d);  // buf_d = ...^(2^64)
        self.mul_unchecked(&buf_d, &mut buf_a);  // buf_a = ... * t1
        buf_a.sqn(30, (&mut buf_b, &mut buf_c), result);  // result = x
    }

    /// Returns self^(2^n) mod p.
    fn sqn(&self, n: usize, scratch: (&mut Self, &mut Self), result: &mut Self) {
        let mut x = scratch.0;
        let mut buffer = scratch.1;

        if n == 1 {
            self.square(result);
            return;
        } else if n == 0 {
            *result = *self;
            return;
        }

        self.square(x);

        let mut i = 2;
        while i < n {
            x.square(buffer);
            i += 1;
            core::mem::swap(&mut x, &mut buffer);
        }

        x.square(result);
    }

    /// Returns self^2 mod p
    pub fn square(&self, result: &mut Self) {
        self.mul_unchecked(self, result);
    }
}

fn bytes_to_u32_words_le_12(bytes: &[u8]) -> [u32; 12] {
    let mut words = [0u32; 12];

    // Process 4 bytes at a time to create little-endian u32 words
    for (i, chunk) in bytes.chunks(4).enumerate() {
        // Convert each big-endian chunk to a little-endian u32
        words[11 - i] = u32::from_be_bytes(chunk.try_into().unwrap());
    }

    words
}

pub fn felt_to_u32_words_le_384<C>(data: &C::FieldElement) -> [u32; 12]
where
    C: PrimeCurveParams384,
{
    bytes_to_u32_words_le_12(data.to_repr().as_slice())
}

#[inline]
fn affine_to_r0_affine_384<C>(affine: &AffinePoint<C>) -> ec::AffinePoint<12, C>
where
    C: PrimeCurveParams384,
{
    if bool::from(affine.is_identity()) {
        return ec::AffinePoint::IDENTITY;
    }

    let x = felt_to_u32_words_le_384::<C>(&affine.x);
    let y = felt_to_u32_words_le_384::<C>(&affine.y);
    ec::AffinePoint::new_unchecked(x, y)
}

pub fn projective_to_affine_384<C>(p: &ProjectivePoint<C>) -> ec::AffinePoint<12, C>
where
    C: PrimeCurveParams384,
{
    let aff = p.to_affine();
    affine_to_r0_affine_384(&aff)
}

pub fn affine_to_projective_384<C>(affine: &ec::AffinePoint<12, C>) -> ProjectivePoint<C>
where
    C: PrimeCurveParams384,
{
    if let Some(value) = affine.as_u32s() {
        // This should only not be within the modulus with a malicious host, panic in that case.
        let x = C::from_u32_words_le(value[0]);
        let y = C::from_u32_words_le(value[1]);

        let affine = AffinePoint { x, y, infinity: 0 };
        ProjectivePoint::from(affine)
    } else {
        ProjectivePoint::IDENTITY
    }
}

pub fn scalar_to_words_384<C>(s: &Scalar<C>) -> [u32; 12]
where
    C: PrimeCurveParams384,
{
    bytes_to_u32_words_le_12(s.to_repr().as_slice())
}

pub mod ec_impl_384 {
    use super::*;

    pub fn mul<C>(lhs: &ProjectivePoint<C>, rhs: &Scalar<C>) -> ProjectivePoint<C>
    where
        C: PrimeCurveParams384,
    {
        let scalar = scalar_to_words_384::<C>(rhs);
        let affine = projective_to_affine_384::<C>(lhs);

        let mut result = ec::AffinePoint::new_unchecked([0u32; 12], [0u32; 12]);
        affine.mul(&scalar, &mut result);
        return affine_to_projective_384(&result);
    }

    pub fn add<C>(lhs: &ProjectivePoint<C>, rhs: &ProjectivePoint<C>) -> ProjectivePoint<C>
    where
        C: PrimeCurveParams384,
    {
        let lhs = projective_to_affine_384::<C>(lhs);
        let rhs = projective_to_affine_384::<C>(rhs);

        let mut result = ec::AffinePoint::new_unchecked([0u32; 12], [0u32; 12]);
        lhs.add(&rhs, &mut result);
        return affine_to_projective_384(&result);
    }

    #[inline]
    pub fn add_mixed<C>(lhs: &ProjectivePoint<C>, rhs: &AffinePoint<C>) -> ProjectivePoint<C>
    where
        C: PrimeCurveParams384,
    {
        let lhs = projective_to_affine_384::<C>(lhs);
        let rhs = affine_to_r0_affine_384(rhs);

        let mut result = ec::AffinePoint::new_unchecked([0u32; 12], [0u32; 12]);
        lhs.add(&rhs, &mut result);
        return affine_to_projective_384(&result);
    }

    pub fn double<C>(point: &ProjectivePoint<C>) -> ProjectivePoint<C>
    where
        C: PrimeCurveParams384,
    {
        let point = projective_to_affine_384::<C>(point);

        let mut result = ec::AffinePoint::new_unchecked([0u32; 12], [0u32; 12]);
        point.double(&mut result);
        return affine_to_projective_384(&result);
    }
}
