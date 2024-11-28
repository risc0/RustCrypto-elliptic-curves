//! A pure-Rust implementation of group operations on secp256k1.

pub(crate) mod affine;
mod field;
#[cfg(feature = "hash2curve")]
mod hash2curve;
mod mul;
pub(crate) mod projective;
pub(crate) mod scalar;

#[cfg(test)]
mod dev;

pub use field::FieldElement;

use self::{affine::AffinePoint, projective::ProjectivePoint, scalar::Scalar};
use crate::Secp256k1;
use elliptic_curve::CurveArithmetic;

impl CurveArithmetic for Secp256k1 {
    type AffinePoint = AffinePoint;
    type ProjectivePoint = ProjectivePoint;
    type Scalar = Scalar;
}

const CURVE_EQUATION_B_SINGLE: u32 = 7u32;

#[rustfmt::skip]
pub(crate) const CURVE_EQUATION_B: FieldElement = FieldElement::from_bytes_unchecked(&[
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, CURVE_EQUATION_B_SINGLE as u8,
]);

#[cfg(all(target_os = "zkvm", target_arch = "riscv32"))]
use risc0_bigint2::ec;

#[cfg(all(target_os = "zkvm", target_arch = "riscv32"))]
use ecdsa_core::elliptic_curve::group::prime::PrimeCurveAffine;

#[cfg(all(target_os = "zkvm", target_arch = "riscv32"))]
pub(crate) fn affine_to_bigint2_affine(
    affine: &AffinePoint,
) -> ec::AffinePoint<8, ec::Secp256k1Curve> {
    if affine.is_identity().into() {
        return ec::AffinePoint::IDENTITY;
    }
    let mut buffer = [[0u32; 8]; 2];
    // TODO this could potentially read from internal repr (check risc0 felt endianness)
    let mut x_bytes: [u8; 32] = affine.x.to_bytes().into();
    let mut y_bytes: [u8; 32] = affine.y.to_bytes().into();
    x_bytes.reverse();
    y_bytes.reverse();

    let x = bytemuck::cast::<_, [u32; 8]>(x_bytes);
    let y = bytemuck::cast::<_, [u32; 8]>(y_bytes);
    ec::AffinePoint::new_unchecked(x, y)
}

#[cfg(all(target_os = "zkvm", target_arch = "riscv32"))]
pub(crate) fn projective_to_affine(p: &ProjectivePoint) -> ec::AffinePoint<8, ec::Secp256k1Curve> {
    let aff = p.to_affine();
    affine_to_bigint2_affine(&aff)
}

#[cfg(all(target_os = "zkvm", target_arch = "riscv32"))]
pub(crate) fn affine_to_projective(
    affine: &ec::AffinePoint<8, ec::Secp256k1Curve>,
) -> ProjectivePoint {
    if let Some(value) = affine.as_u32s() {
        let mut x = bytemuck::cast::<_, [u8; 32]>(value[0]);
        let mut y = bytemuck::cast::<_, [u8; 32]>(value[1]);
        x.reverse();
        y.reverse();

        crate::AffinePoint::new(
            FieldElement::from_bytes_unchecked(&x),
            FieldElement::from_bytes_unchecked(&y),
        )
        .into()
    } else {
        ProjectivePoint::IDENTITY
    }
}

#[cfg(all(target_os = "zkvm", target_arch = "riscv32"))]
pub(crate) fn scalar_to_words(s: &Scalar) -> [u32; 8] {
    let mut bytes: [u8; 32] = s.to_bytes().into();
    // U256 is big endian, need to flip to little endian.
    bytes.reverse();
    bytemuck::cast::<_, [u32; 8]>(bytes)
}

#[cfg(test)]
mod tests {
    use super::CURVE_EQUATION_B;
    use hex_literal::hex;

    const CURVE_EQUATION_B_BYTES: [u8; 32] =
        hex!("0000000000000000000000000000000000000000000000000000000000000007");

    #[test]
    fn verify_constants() {
        assert_eq!(CURVE_EQUATION_B.to_bytes(), CURVE_EQUATION_B_BYTES.into());
    }

    #[test]
    fn generate_secret_key() {
        use crate::SecretKey;
        use elliptic_curve::rand_core::OsRng;
        let key = SecretKey::random(&mut OsRng);

        // Sanity check
        assert!(!key.to_bytes().iter().all(|b| *b == 0))
    }
}
