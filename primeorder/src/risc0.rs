use crate::FieldBytes;
use core::marker::PhantomData;
use core::ops::Deref;
use elliptic_curve::{bigint::ArrayEncoding, generic_array::GenericArray, PrimeField};

use crate::PrimeCurveParams;

/// Representation of a field element in raw bytes form. This is not in montgomery form.
#[derive(Copy, Clone, Default, Debug)]
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
    pub fn mul_unchecked(&self, rhs: &Self, result: &mut Self) {
        risc0_bigint2::field::modmul_256_unchecked(
            &self.data,
            &rhs.data,
            &C::PRIME_LE_WORDS,
            &mut result.data,
        );
    }

    pub fn add_unchecked(&self, rhs: &Self, result: &mut Self) {
        risc0_bigint2::field::modadd_256_unchecked(
            &self.data,
            &rhs.data,
            &C::PRIME_LE_WORDS,
            &mut result.data,
        );
    }
}
