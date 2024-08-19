//! Pure Rust implementation of group operations on secp256r1.
//!
//! Curve parameters can be found in [NIST SP 800-186] § G.1.2: Curve P-256.
//!
//! [NIST SP 800-186]: https://csrc.nist.gov/publications/detail/sp/800-186/final

pub(crate) mod affine;
pub(crate) mod field;
#[cfg(feature = "hash2curve")]
mod hash2curve;
pub(crate) mod projective;
pub(crate) mod scalar;
pub(crate) mod util;

use self::{field::FieldElement, scalar::Scalar};
use crate::NistP256;
use elliptic_curve::{bigint::U256, CurveArithmetic};

/// Elliptic curve point in affine coordinates.
pub type AffinePoint = affine::AffinePoint;

/// Elliptic curve point in projective coordinates.
pub type ProjectivePoint = projective::ProjectivePoint;

impl CurveArithmetic for NistP256 {
    type AffinePoint = AffinePoint;
    type ProjectivePoint = ProjectivePoint;
    type Scalar = Scalar;
}

/// a = -3
const CURVE_EQUATION_A: FieldElement = FieldElement(U256::from_be_hex(
    "FFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFC",
));

const CURVE_EQUATION_B: FieldElement = FieldElement(U256::from_be_hex(
    "5AC635D8AA3A93E7B3EBBD55769886BC651D06B0CC53B0F63BCE3C3E27D2604B",
));

#[cfg(test)]
mod tests {
    use super::FieldElement;
    use crate::{
        arithmetic::{CURVE_EQUATION_A, CURVE_EQUATION_B},
        AffinePoint,
    };

    #[test]
    fn equation_a_constant() {
        let equation_a = FieldElement::from_u64(3).neg();
        assert_eq!(equation_a, CURVE_EQUATION_A);
    }

    #[test]
    fn equation_b_constant() {
        let equation_b = FieldElement::from_hex(
            "5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b",
        );
        assert_eq!(equation_b, CURVE_EQUATION_B);
    }

    #[test]
    fn generator_constant() {
        let generator: (FieldElement, FieldElement) = (
            FieldElement::from_hex(
                "6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296",
            ),
            FieldElement::from_hex(
                "4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5",
            ),
        );
        assert_eq!(generator.0, AffinePoint::GENERATOR.x);
        assert_eq!(generator.1, AffinePoint::GENERATOR.y);
    }
}
