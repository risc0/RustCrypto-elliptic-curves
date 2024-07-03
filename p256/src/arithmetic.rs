//! Pure Rust implementation of group operations on secp256r1.
//!
//! Curve parameters can be found in [NIST SP 800-186] § G.1.2: Curve P-256.
//!
//! [NIST SP 800-186]: https://csrc.nist.gov/publications/detail/sp/800-186/final

pub(crate) mod field;
#[cfg(feature = "hash2curve")]
mod hash2curve;
pub(crate) mod scalar;

use self::{field::FieldElement, scalar::Scalar};
use crate::NistP256;
use elliptic_curve::{bigint::U256, CurveArithmetic, PrimeCurveArithmetic};
use primeorder::{point_arithmetic, PrimeCurveParams};

/// Elliptic curve point in affine coordinates.
pub type AffinePoint = primeorder::AffinePoint<NistP256>;

/// Elliptic curve point in projective coordinates.
pub type ProjectivePoint = primeorder::ProjectivePoint<NistP256>;

impl CurveArithmetic for NistP256 {
    type AffinePoint = AffinePoint;
    type ProjectivePoint = ProjectivePoint;
    type Scalar = Scalar;
}

impl PrimeCurveArithmetic for NistP256 {
    type CurveGroup = ProjectivePoint;
}

/// Adapted from [NIST SP 800-186] § G.1.2: Curve P-256.
///
/// [NIST SP 800-186]: https://csrc.nist.gov/publications/detail/sp/800-186/final
impl PrimeCurveParams for NistP256 {
    type FieldElement = FieldElement;
    type PointArithmetic = point_arithmetic::EquationAIsMinusThree;

    /// a = -3
    const EQUATION_A: FieldElement = FieldElement(U256::from_be_hex(
        "FFFFFFFC00000004000000000000000000000003FFFFFFFFFFFFFFFFFFFFFFFC",
    ));

    const EQUATION_B: FieldElement = FieldElement(U256::from_be_hex(
        "DC30061D04874834E5A220ABF7212ED6ACF005CD78843090D89CDF6229C4BDDF",
    ));

    /// Base point of P-256.
    ///
    /// Defined in NIST SP 800-186 § G.1.2:
    ///
    /// ```text
    /// Gₓ = 6b17d1f2 e12c4247 f8bce6e5 63a440f2 77037d81 2deb33a0 f4a13945 d898c296
    /// Gᵧ = 4fe342e2 fe1a7f9b 8ee7eb4a 7c0f9e16 2bce3357 6b315ece cbb64068 37bf51f5
    /// ```
    const GENERATOR: (FieldElement, FieldElement) = (
        FieldElement(U256::from_be_hex(
            "18905F76A53755C679FB732B7762251075BA95FC5FEDB60179E730D418A9143C",
        )),
        FieldElement(U256::from_be_hex(
            "8571FF1825885D85D2E88688DD21F3258B4AB8E4BA19E45CDDF25357CE95560A",
        )),
    );
}

#[cfg(test)]
mod tests {
    use super::FieldElement;
    use crate::NistP256;
    use primeorder::PrimeCurveParams;

    #[test]
    fn equation_a_constant() {
        let equation_a = FieldElement::from_u64(3).neg();
        assert_eq!(equation_a, NistP256::EQUATION_A);
    }

    #[test]
    fn equation_b_constant() {
        let equation_b = FieldElement::from_hex(
            "5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b",
        );
        assert_eq!(equation_b, NistP256::EQUATION_B);
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
        assert_eq!(generator, NistP256::GENERATOR);
    }
}
