//! Affine curve points.

#![allow(clippy::op_ref)]

use super::{FieldElement, ProjectivePoint, CURVE_EQUATION_A, CURVE_EQUATION_B};
use crate::{CompressedPoint, EncodedPoint, FieldBytes, NistP256, PublicKey, Scalar};
use core::ops::{Mul, Neg};
use elliptic_curve::{
    bigint::U256,
    ff::PrimeField,
    group::{prime::PrimeCurveAffine, GroupEncoding},
    point::{AffineCoordinates, DecompactPoint, DecompressPoint},
    sec1::{self, FromEncodedPoint, ToCompactEncodedPoint, ToEncodedPoint},
    subtle::{Choice, ConditionallySelectable, ConstantTimeEq, ConstantTimeGreater, CtOption},
    zeroize::DefaultIsZeroes,
    Error, FieldBytesEncoding, Result,
};

#[cfg(feature = "serde")]
use serdect::serde::{de, ser, Deserialize, Serialize};

/// Point on a Weierstrass curve in affine coordinates.
#[derive(Clone, Copy, Debug)]
pub struct AffinePoint {
    /// x-coordinate
    pub(crate) x: FieldElement,

    /// y-coordinate
    pub(crate) y: FieldElement,

    /// Is this point the point at infinity? 0 = no, 1 = yes
    ///
    /// This is a proxy for [`Choice`], but uses `u8` instead to permit `const`
    /// constructors for `IDENTITY` and `GENERATOR`.
    pub(crate) infinity: u8,
}

// NOTE: Removes the PrimeCurveParams generic. TODO impl PrimeCurveParams
impl AffinePoint {
    /// Additive identity of the group a.k.a. the point at infinity.
    pub const IDENTITY: Self = Self {
        x: FieldElement::ZERO,
        y: FieldElement::ZERO,
        infinity: 1,
    };

    /// Base point of the curve.
    pub const GENERATOR: Self = Self {
        x: FieldElement(U256::from_be_hex(
            "6B17D1F2E12C4247F8BCE6E563A440F277037D812DEB33A0F4A13945D898C296",
        )),
        y: FieldElement(U256::from_be_hex(
            "4FE342E2FE1A7F9B8EE7EB4A7C0F9E162BCE33576B315ECECBB6406837BF51F5",
        )),
        infinity: 0,
    };

    /// Is this point the point at infinity?
    pub fn is_identity(&self) -> Choice {
        Choice::from(self.infinity)
    }

    /// Conditionally negate [`AffinePoint`] for use with point compaction.
    fn to_compact(self) -> Self {
        let neg_self = -self;
        let choice = <U256 as FieldBytesEncoding<NistP256>>::decode_field_bytes(&self.y.to_repr())
            .ct_gt(&<U256 as FieldBytesEncoding<NistP256>>::decode_field_bytes(
                &neg_self.y.to_repr(),
            ));

        Self {
            x: self.x,
            y: FieldElement::conditional_select(&self.y, &neg_self.y, choice),
            infinity: self.infinity,
        }
    }
}

impl AffineCoordinates for AffinePoint {
    type FieldRepr = FieldBytes;

    fn x(&self) -> FieldBytes {
        self.x.to_repr()
    }

    fn y_is_odd(&self) -> Choice {
        self.y.is_odd()
    }
}

impl ConditionallySelectable for AffinePoint {
    #[inline(always)]
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        Self {
            x: FieldElement::conditional_select(&a.x, &b.x, choice),
            y: FieldElement::conditional_select(&a.y, &b.y, choice),
            infinity: u8::conditional_select(&a.infinity, &b.infinity, choice),
        }
    }
}

impl ConstantTimeEq for AffinePoint {
    fn ct_eq(&self, other: &Self) -> Choice {
        self.x.ct_eq(&other.x) & self.y.ct_eq(&other.y) & self.infinity.ct_eq(&other.infinity)
    }
}

impl Default for AffinePoint {
    fn default() -> Self {
        Self::IDENTITY
    }
}

impl DefaultIsZeroes for AffinePoint {}

impl DecompressPoint<NistP256> for AffinePoint {
    fn decompress(x_bytes: &FieldBytes, y_is_odd: Choice) -> CtOption<Self> {
        FieldElement::from_repr(*x_bytes).and_then(|x| {
            let alpha = x * &x * &x + &(CURVE_EQUATION_A * &x) + &CURVE_EQUATION_B;
            let beta = alpha.sqrt();

            beta.map(|beta| {
                let y =
                    FieldElement::conditional_select(&-beta, &beta, beta.is_odd().ct_eq(&y_is_odd));

                Self { x, y, infinity: 0 }
            })
        })
    }
}

impl DecompactPoint<NistP256> for AffinePoint {
    fn decompact(x_bytes: &FieldBytes) -> CtOption<Self> {
        Self::decompress(x_bytes, Choice::from(0)).map(|point| point.to_compact())
    }
}

impl FromEncodedPoint<NistP256> for AffinePoint {
    /// Attempts to parse the given [`EncodedPoint`] as an SEC1-encoded
    /// [`AffinePoint`].
    ///
    /// # Returns
    ///
    /// `None` value if `encoded_point` is not on the secp384r1 curve.
    fn from_encoded_point(encoded_point: &EncodedPoint) -> CtOption<Self> {
        match encoded_point.coordinates() {
            sec1::Coordinates::Identity => CtOption::new(Self::IDENTITY, 1.into()),
            sec1::Coordinates::Compact { x } => Self::decompact(x),
            sec1::Coordinates::Compressed { x, y_is_odd } => {
                Self::decompress(x, Choice::from(y_is_odd as u8))
            }
            sec1::Coordinates::Uncompressed { x, y } => FieldElement::from_repr(*y).and_then(|y| {
                FieldElement::from_repr(*x).and_then(|x| {
                    let lhs = y * &y;
                    let rhs = x * &x * &x + &(CURVE_EQUATION_A * &x) + &CURVE_EQUATION_B;
                    CtOption::new(Self { x, y, infinity: 0 }, lhs.ct_eq(&rhs))
                })
            }),
        }
    }
}

impl Eq for AffinePoint {}

impl From<ProjectivePoint> for AffinePoint {
    fn from(p: ProjectivePoint) -> AffinePoint {
        p.to_affine()
    }
}

impl From<&ProjectivePoint> for AffinePoint {
    fn from(p: &ProjectivePoint) -> AffinePoint {
        p.to_affine()
    }
}

impl From<PublicKey> for AffinePoint {
    fn from(public_key: PublicKey) -> AffinePoint {
        *public_key.as_affine()
    }
}

impl From<&PublicKey> for AffinePoint {
    fn from(public_key: &PublicKey) -> AffinePoint {
        AffinePoint::from(*public_key)
    }
}

impl From<AffinePoint> for EncodedPoint {
    fn from(affine: AffinePoint) -> EncodedPoint {
        affine.to_encoded_point(false)
    }
}

impl GroupEncoding for AffinePoint {
    type Repr = CompressedPoint;

    /// NOTE: not constant-time with respect to identity point
    fn from_bytes(bytes: &Self::Repr) -> CtOption<Self> {
        EncodedPoint::from_bytes(bytes)
            .map(|point| CtOption::new(point, Choice::from(1)))
            .unwrap_or_else(|_| {
                // SEC1 identity encoding is technically 1-byte 0x00, but the
                // `GroupEncoding` API requires a fixed-width `Repr`
                let is_identity = bytes.ct_eq(&Self::Repr::default());
                CtOption::new(EncodedPoint::identity(), is_identity)
            })
            .and_then(|point| Self::from_encoded_point(&point))
    }

    fn from_bytes_unchecked(bytes: &Self::Repr) -> CtOption<Self> {
        // No unchecked conversion possible for compressed points
        Self::from_bytes(bytes)
    }

    fn to_bytes(&self) -> Self::Repr {
        let encoded = self.to_encoded_point(true);
        let mut result = CompressedPoint::default();
        result[..encoded.len()].copy_from_slice(encoded.as_bytes());
        result
    }
}

impl PartialEq for AffinePoint {
    fn eq(&self, other: &Self) -> bool {
        self.ct_eq(other).into()
    }
}

impl PrimeCurveAffine for AffinePoint {
    type Curve = ProjectivePoint;
    type Scalar = Scalar;

    fn identity() -> AffinePoint {
        Self::IDENTITY
    }

    fn generator() -> AffinePoint {
        Self::GENERATOR
    }

    fn is_identity(&self) -> Choice {
        self.is_identity()
    }

    fn to_curve(&self) -> ProjectivePoint {
        ProjectivePoint::from(*self)
    }
}

impl ToCompactEncodedPoint<NistP256> for AffinePoint {
    /// Serialize this value as a  SEC1 compact [`EncodedPoint`]
    fn to_compact_encoded_point(&self) -> CtOption<EncodedPoint> {
        let point = self.to_compact();

        let mut bytes = CompressedPoint::default();
        bytes[0] = sec1::Tag::Compact.into();
        bytes[1..].copy_from_slice(&point.x.to_repr());

        let encoded = EncodedPoint::from_bytes(bytes);
        let is_some = point.y.ct_eq(&self.y);
        CtOption::new(encoded.unwrap_or_default(), is_some)
    }
}

impl ToEncodedPoint<NistP256> for AffinePoint {
    fn to_encoded_point(&self, compress: bool) -> EncodedPoint {
        EncodedPoint::conditional_select(
            &EncodedPoint::from_affine_coordinates(&self.x.to_repr(), &self.y.to_repr(), compress),
            &EncodedPoint::identity(),
            self.is_identity(),
        )
    }
}

impl TryFrom<EncodedPoint> for AffinePoint {
    type Error = Error;

    fn try_from(point: EncodedPoint) -> Result<AffinePoint> {
        AffinePoint::try_from(&point)
    }
}

impl TryFrom<&EncodedPoint> for AffinePoint {
    type Error = Error;

    fn try_from(point: &EncodedPoint) -> Result<AffinePoint> {
        Option::from(AffinePoint::from_encoded_point(point)).ok_or(Error)
    }
}

impl TryFrom<AffinePoint> for PublicKey {
    type Error = Error;

    fn try_from(affine_point: AffinePoint) -> Result<PublicKey> {
        PublicKey::from_affine(affine_point)
    }
}

impl TryFrom<&AffinePoint> for PublicKey {
    type Error = Error;

    fn try_from(affine_point: &AffinePoint) -> Result<PublicKey> {
        PublicKey::try_from(*affine_point)
    }
}

//
// Arithmetic trait impls
//

impl Mul<Scalar> for AffinePoint {
    type Output = ProjectivePoint;

    fn mul(self, scalar: Scalar) -> ProjectivePoint {
        ProjectivePoint::from(self) * scalar
    }
}

impl Mul<&Scalar> for AffinePoint {
    type Output = ProjectivePoint;

    fn mul(self, scalar: &Scalar) -> ProjectivePoint {
        ProjectivePoint::from(self) * *scalar
    }
}

impl Neg for AffinePoint {
    type Output = Self;

    fn neg(self) -> Self {
        AffinePoint {
            x: self.x,
            y: -self.y,
            infinity: self.infinity,
        }
    }
}

impl Neg for &AffinePoint {
    type Output = AffinePoint;

    fn neg(self) -> AffinePoint {
        -(*self)
    }
}

//
// serde support
//

#[cfg(feature = "serde")]
impl Serialize for AffinePoint {
    fn serialize<S>(&self, serializer: S) -> core::result::Result<S::Ok, S::Error>
    where
        S: ser::Serializer,
    {
        self.to_encoded_point(true).serialize(serializer)
    }
}

#[cfg(feature = "serde")]
impl<'de> Deserialize<'de> for AffinePoint {
    fn deserialize<D>(deserializer: D) -> core::result::Result<Self, D::Error>
    where
        D: de::Deserializer<'de>,
    {
        EncodedPoint::deserialize(deserializer)?
            .try_into()
            .map_err(de::Error::custom)
    }
}
