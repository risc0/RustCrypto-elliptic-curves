//! Pure Rust implementation of group operations on secp256r1.
//!
//! Curve parameters can be found in [NIST SP 800-186] § G.1.2: Curve P-256.
//!
//! [NIST SP 800-186]: https://csrc.nist.gov/publications/detail/sp/800-186/final

pub(crate) mod field;
#[cfg(feature = "hash2curve")]
mod hash2curve;
pub(crate) mod scalar;
pub(crate) mod util;

use self::{field::FieldElement, scalar::Scalar};
use crate::NistP256;
use elliptic_curve::{CurveArithmetic, PrimeCurveArithmetic};
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

#[cfg(all(target_os = "zkvm", target_arch = "riscv32"))]
use primeorder::__risc0::FieldElement256;

/// Adapted from [NIST SP 800-186] § G.1.2: Curve P-256.
///
/// [NIST SP 800-186]: https://csrc.nist.gov/publications/detail/sp/800-186/final
impl PrimeCurveParams for NistP256 {
    type FieldElement = FieldElement;
    type PointArithmetic = point_arithmetic::EquationAIsMinusThree;

    /// a = -3
    const EQUATION_A: FieldElement = FieldElement::from_u64(3).neg();

    const EQUATION_B: FieldElement =
        FieldElement::from_hex("5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b");

    /// Base point of P-256.
    ///
    /// Defined in NIST SP 800-186 § G.1.2:
    ///
    /// ```text
    /// Gₓ = 6b17d1f2 e12c4247 f8bce6e5 63a440f2 77037d81 2deb33a0 f4a13945 d898c296
    /// Gᵧ = 4fe342e2 fe1a7f9b 8ee7eb4a 7c0f9e16 2bce3357 6b315ece cbb64068 37bf51f5
    /// ```
    const GENERATOR: (FieldElement, FieldElement) = (
        FieldElement::from_hex("6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296"),
        FieldElement::from_hex("4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5"),
    );

    #[cfg(all(target_os = "zkvm", target_arch = "riscv32"))]
    const PRIME_LE_WORDS: [u32; 8] = crate::__risc0::SECP256R1_PRIME;

    #[cfg(all(target_os = "zkvm", target_arch = "riscv32"))]
    const ORDER_LE_WORDS: [u32; 8] = crate::__risc0::SECP256R1_ORDER;

    #[cfg(all(target_os = "zkvm", target_arch = "riscv32"))]
    const EQUATION_A_LE: FieldElement256<NistP256> =
        FieldElement256::new_unchecked(crate::__risc0::SECP256R1_EQUATION_A_LE);

    #[cfg(all(target_os = "zkvm", target_arch = "riscv32"))]
    const EQUATION_B_LE: FieldElement256<NistP256> =
        FieldElement256::new_unchecked(crate::__risc0::SECP256R1_EQUATION_B_LE);

    #[cfg(all(target_os = "zkvm", target_arch = "riscv32"))]
    fn from_u32_words_le(words: [u32; 8]) -> elliptic_curve::subtle::CtOption<FieldElement> {
        FieldElement::from_words_le(words)
    }
}
