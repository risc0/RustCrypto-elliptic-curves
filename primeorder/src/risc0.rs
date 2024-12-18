use crate::FieldBytes;
use crate::{affine::AffinePoint, projective::ProjectivePoint};
use core::marker::PhantomData;
use core::ops::Deref;
use elliptic_curve::generic_array::GenericArray;
use elliptic_curve::subtle::{Choice, ConditionallySelectable};
use elliptic_curve::{PrimeField, Scalar};
use risc0_bigint2::ec;

use crate::PrimeCurveParams;

/// Representation of a field element in raw bytes form. This is not in montgomery form.
#[derive(Copy, Clone, Default, Debug, PartialEq, Eq)]
pub struct FieldElement256<C> {
    pub(crate) data: [u32; 8],
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
    C: PrimeCurveParams,
{
    pub(crate) fn mul_unchecked(&self, rhs: &Self, result: &mut Self) {
        risc0_bigint2::field::modmul_256_unchecked(
            &self.data,
            &rhs.data,
            &C::PRIME_LE_WORDS,
            &mut result.data,
        );
    }

    pub(crate) fn add_unchecked(&self, rhs: &Self, result: &mut Self) {
        risc0_bigint2::field::modadd_256_unchecked(
            &self.data,
            &rhs.data,
            &C::PRIME_LE_WORDS,
            &mut result.data,
        );
    }

    /// Calculate the square root of the field element, using the provided buffer as scratch, and
    /// writing the result to the `result` parameter.
    pub(crate) fn sqrt_unchecked(&self, scratch: &mut Self, result: &mut Self) {
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
    pub(crate) fn square(&self, result: &mut Self) {
        self.mul_unchecked(self, result);
    }
}

// TODO remove inline
#[inline(never)]
pub(crate) fn projective_to_affine<C>(p: &ProjectivePoint<C>) -> ec::AffinePoint<8, C>
where
    C: PrimeCurveParams,
{
    let aff = p.to_affine();
    let x_bytes = aff.x.to_repr();
    let y_bytes = aff.y.to_repr();
    let mut x_bytes_arr: [u8; 32] = x_bytes.as_slice().try_into().unwrap();
    let mut y_bytes_arr: [u8; 32] = y_bytes.as_slice().try_into().unwrap();
    x_bytes_arr.reverse();
    y_bytes_arr.reverse();
    let x = bytemuck::cast::<_, [u32; 8]>(x_bytes_arr);
    let y = bytemuck::cast::<_, [u32; 8]>(y_bytes_arr);
    ec::AffinePoint::new_unchecked(x, y)
}

// TODO remove inline
#[inline(never)]
pub(crate) fn affine_to_projective<C>(affine: &ec::AffinePoint<8, C>) -> ProjectivePoint<C>
where
    C: PrimeCurveParams,
{
    if let Some(value) = affine.as_u32s() {
        // TODO a lot of potentially unnecessary copying here.
        let mut x = bytemuck::cast::<_, [u8; 32]>(value[0]);
        let mut y = bytemuck::cast::<_, [u8; 32]>(value[1]);
        x.reverse();
        y.reverse();
        let x_arr = GenericArray::from_slice(&x);
        let y_arr = GenericArray::from_slice(&y);
        let affine = AffinePoint {
            x: C::FieldElement::from_repr(x_arr.clone()).unwrap(),
            y: C::FieldElement::from_repr(y_arr.clone()).unwrap(),
            infinity: 0,
        };
        ProjectivePoint::from(affine)
    } else {
        ProjectivePoint::IDENTITY
    }
}

#[inline(never)]
pub(crate) fn scalar_to_words<C>(s: &Scalar<C>) -> [u32; 8]
where
    C: PrimeCurveParams,
{
    let mut bytes: [u8; 32] = s.to_repr().as_slice().try_into().unwrap();
    // U256 is big endian, need to flip to little endian.
    bytes.reverse();
    bytemuck::cast::<_, [u32; 8]>(bytes)
}

pub(crate) mod ec_impl {
    use super::*;

    pub(crate) fn mul<C>(lhs: &ProjectivePoint<C>, rhs: &Scalar<C>) -> ProjectivePoint<C>
    where
        C: PrimeCurveParams,
    {
        let scalar = scalar_to_words::<C>(rhs);
        let affine = projective_to_affine::<C>(lhs);

        let mut result = risc0_bigint2::ec::AffinePoint::new_unchecked([0u32; 8], [0u32; 8]);
        affine.mul(&scalar, &mut result);
        return affine_to_projective(&result);
    }

    pub(crate) fn add<C>(lhs: &ProjectivePoint<C>, rhs: &ProjectivePoint<C>) -> ProjectivePoint<C>
    where
        C: PrimeCurveParams,
    {
        let lhs = projective_to_affine::<C>(lhs);
        let rhs = projective_to_affine::<C>(rhs);

        let mut result = risc0_bigint2::ec::AffinePoint::new_unchecked([0u32; 8], [0u32; 8]);
        lhs.add(&rhs, &mut result);
        return affine_to_projective(&result);
    }

    /// Implements complete mixed addition for curves with `a = -3`
    ///
    /// Implements the complete mixed addition formula from [Renes-Costello-Batina 2015]
    /// (Algorithm 5). The comments after each line indicate which algorithm
    /// steps are being performed.
    ///
    /// [Renes-Costello-Batina 2015]: https://eprint.iacr.org/2015/1060
    pub(crate) fn add_mixed<C>(lhs: &ProjectivePoint<C>, rhs: &AffinePoint<C>) -> ProjectivePoint<C>
    where
        C: PrimeCurveParams,
    {
        // TODO
        todo!()
    }

    /// Implements point doubling for curves with `a = -3`
    ///
    /// Implements the exception-free point doubling formula from [Renes-Costello-Batina 2015]
    /// (Algorithm 6). The comments after each line indicate which algorithm
    /// steps are being performed.
    ///
    /// [Renes-Costello-Batina 2015]: https://eprint.iacr.org/2015/1060
    pub(crate) fn double<C>(point: &ProjectivePoint<C>) -> ProjectivePoint<C>
    where
        C: PrimeCurveParams,
    {
        // TODO
        todo!()
    }
}
