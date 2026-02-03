#![no_std]
#![cfg_attr(docsrs, feature(doc_auto_cfg))]
#![doc(
    html_logo_url = "https://raw.githubusercontent.com/RustCrypto/meta/master/logo.svg",
    html_favicon_url = "https://raw.githubusercontent.com/RustCrypto/meta/master/logo.svg"
)]
#![forbid(unsafe_code)]
#![warn(missing_docs, rust_2018_idioms, unused_qualifications)]
#![doc = include_str!("../README.md")]

pub mod point_arithmetic;

mod affine;
#[cfg(feature = "dev")]
mod dev;
mod field;
mod projective;

/// WARNING: This is not part of the public API of the crate, and will not follow semver guarantees.
#[cfg(all(target_os = "zkvm", target_arch = "riscv32"))]
#[path = "risc0.rs"]
pub mod __risc0;

pub use crate::{affine::AffinePoint, projective::ProjectivePoint};
pub use elliptic_curve::{self, point::Double, Field, FieldBytes, PrimeCurve, PrimeField};

use elliptic_curve::CurveArithmetic;

/// Parameters for elliptic curves of prime order which can be described by the
/// short Weierstrass equation.
#[cfg(not(all(target_os = "zkvm", target_arch = "riscv32")))]
pub trait PrimeCurveParams:
    PrimeCurve
    + CurveArithmetic
    + CurveArithmetic<AffinePoint = AffinePoint<Self>>
    + CurveArithmetic<ProjectivePoint = ProjectivePoint<Self>>
{
    /// Base field element type.
    type FieldElement: PrimeField<Repr = FieldBytes<Self>>;

    /// [Point arithmetic](point_arithmetic) implementation, might be optimized for this specific curve
    type PointArithmetic: point_arithmetic::PointArithmetic<Self>;

    /// Coefficient `a` in the curve equation.
    const EQUATION_A: Self::FieldElement;

    /// Coefficient `b` in the curve equation.
    const EQUATION_B: Self::FieldElement;

    /// Generator point's affine coordinates: (x, y).
    const GENERATOR: (Self::FieldElement, Self::FieldElement);
}

/// Parameters for 384-bit elliptic curves (non-zkvm version).
/// Alias for PrimeCurveParams since they're identical on non-zkvm.
#[cfg(not(all(target_os = "zkvm", target_arch = "riscv32")))]
pub trait PrimeCurveParams384: PrimeCurveParams {}

/// Blanket implementation: any PrimeCurveParams is also PrimeCurveParams384 on non-zkvm.
#[cfg(not(all(target_os = "zkvm", target_arch = "riscv32")))]
impl<T: PrimeCurveParams> PrimeCurveParams384 for T {}

#[cfg(all(target_os = "zkvm", target_arch = "riscv32"))]
use risc0_bigint2::ec;

/// Parameters for elliptic curves of prime order which can be described by the
/// short Weierstrass equation.
#[cfg(all(target_os = "zkvm", target_arch = "riscv32"))]
pub trait PrimeCurveParams:
    Sized +
    PrimeCurve
    + CurveArithmetic
    + CurveArithmetic<AffinePoint = AffinePoint<Self>>
    + CurveArithmetic<ProjectivePoint = ProjectivePoint<Self>>
{
    /// Base field element type.
    type FieldElement: PrimeField<Repr = FieldBytes<Self>>;

    /// [Point arithmetic](point_arithmetic) implementation, might be optimized for this specific curve
    type PointArithmetic: point_arithmetic::PointArithmetic<Self>;

    /// Coefficient `a` in the curve equation.
    const EQUATION_A: Self::FieldElement;

    /// Coefficient `b` in the curve equation.
    const EQUATION_B: Self::FieldElement;

    /// Generator point's affine coordinates: (x, y).
    const GENERATOR: (Self::FieldElement, Self::FieldElement);

    /// Scalar multiplication for zkvm (accelerated). Curves implement this based on their width.
    fn zkvm_mul(point: &ProjectivePoint<Self>, scalar: &elliptic_curve::Scalar<Self>) -> ProjectivePoint<Self>;

    /// Convert projective to affine for zkvm (accelerated). Curves implement this based on their width.
    fn zkvm_to_affine(point: &ProjectivePoint<Self>) -> AffinePoint<Self>;

    /// Decompress point for zkvm (accelerated). Curves implement this based on their width.
    fn zkvm_decompress(x_bytes: &FieldBytes<Self>, y_is_odd: elliptic_curve::subtle::Choice) -> elliptic_curve::subtle::CtOption<AffinePoint<Self>>;
}

/// Parameters for 256-bit elliptic curves on zkvm.
#[cfg(all(target_os = "zkvm", target_arch = "riscv32"))]
pub trait PrimeCurveParams256:
    PrimeCurveParams +
    ec::Curve<{ec::EC_256_WIDTH_WORDS}>
{
    /// [Point arithmetic](point_arithmetic) implementation
    type PointArithmetic: point_arithmetic::PointArithmetic<Self>;

    /// Curve prime in little-endian words to be compatible with risc0 expected layout.
    const PRIME_LE_WORDS: [u32; 8];

    /// Order of the curve in little-endian words to be compatible with risc0 expected layout.
    const ORDER_LE_WORDS: [u32; 8];

    /// Coefficient `a` in the curve equation in little-endian words to be compatible with risc0
    /// expected layout.
    const EQUATION_A_LE: __risc0::FieldElement256<Self>;

    /// Coefficient `b` in the curve equation in little-endian words to be compatible with risc0
    /// expected layout.
    const EQUATION_B_LE: __risc0::FieldElement256<Self>;

    /// Convert from little-endian words to field element.
    fn from_u32_words_le(words: [u32; 8]) -> Self::FieldElement;
}

/// Parameters for 384-bit elliptic curves on zkvm.
#[cfg(all(target_os = "zkvm", target_arch = "riscv32"))]
pub trait PrimeCurveParams384:
    PrimeCurveParams +
    ec::Curve<{ec::EC_384_WIDTH_WORDS}>
{
    /// [Point arithmetic](point_arithmetic) implementation
    type PointArithmetic: point_arithmetic::PointArithmetic384<Self>;

    /// Curve prime in little-endian words to be compatible with risc0 expected layout.
    const PRIME_LE_WORDS: [u32; 12];

    /// Order of the curve in little-endian words to be compatible with risc0 expected layout.
    const ORDER_LE_WORDS: [u32; 12];

    /// Coefficient `a` in the curve equation in little-endian words to be compatible with risc0
    /// expected layout.
    const EQUATION_A_LE: __risc0::FieldElement384<Self>;

    /// Coefficient `b` in the curve equation in little-endian words to be compatible with risc0
    /// expected layout.
    const EQUATION_B_LE: __risc0::FieldElement384<Self>;

    /// Convert from little-endian words to field element.
    fn from_u32_words_le(words: [u32; 12]) -> Self::FieldElement;
}
