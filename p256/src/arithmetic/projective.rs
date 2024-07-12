//! Projective curve points.

#![allow(clippy::needless_range_loop, clippy::op_ref)]

use super::{AffinePoint, FieldElement};
use crate::{
    arithmetic::{CURVE_EQUATION_A, CURVE_EQUATION_B},
    CompressedPoint, EncodedPoint, NistP256, PublicKey, Scalar,
};
use core::{
    borrow::Borrow,
    iter::Sum,
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};
use elliptic_curve::{
    bigint::{ArrayEncoding, U256},
    generic_array::ArrayLength,
    group::{
        self,
        cofactor::CofactorGroup,
        ff::Field,
        prime::{PrimeCurve, PrimeGroup},
        Group, GroupEncoding,
    },
    ops::{BatchInvert, LinearCombination, MulByGenerator},
    point::Double,
    rand_core::RngCore,
    sec1::{FromEncodedPoint, ModulusSize, ToEncodedPoint, UncompressedPointSize},
    subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption},
    zeroize::DefaultIsZeroes,
    BatchNormalize, Error, FieldBytes, FieldBytesSize, Result,
};

#[cfg(feature = "alloc")]
use alloc::vec::Vec;

/// Point on a Weierstrass curve in projective coordinates.
#[derive(Clone, Copy, Debug)]
pub struct ProjectivePoint {
    pub(crate) x: FieldElement,
    pub(crate) y: FieldElement,
    pub(crate) z: FieldElement,
}

impl ProjectivePoint {
    /// Additive identity of the group a.k.a. the point at infinity.
    pub const IDENTITY: Self = Self {
        x: FieldElement::ZERO,
        y: FieldElement::ONE,
        z: FieldElement::ZERO,
    };

    /// Base point of the curve.
    pub const GENERATOR: Self = Self {
        x: AffinePoint::GENERATOR.x,
        y: AffinePoint::GENERATOR.y,
        z: FieldElement::ONE,
    };

    /// Returns the affine representation of this point, or `None` if it is the identity.
    pub fn to_affine(&self) -> AffinePoint {
        <FieldElement as Field>::invert(&self.z)
            .map(|zinv| self.to_affine_internal(zinv))
            .unwrap_or(AffinePoint::IDENTITY)
    }

    pub(super) fn to_affine_internal(self, zinv: FieldElement) -> AffinePoint {
        AffinePoint {
            x: self.x * &zinv,
            y: self.y * &zinv,
            infinity: 0,
        }
    }

    /// Returns `-self`.
    pub fn neg(&self) -> Self {
        Self {
            x: self.x,
            y: -self.y,
            z: self.z,
        }
    }

    /// Returns `self + other`.
    /// Implements complete addition for curves with `a = -3`
    ///
    /// Implements the complete addition formula from [Renes-Costello-Batina 2015]
    /// (Algorithm 4). The comments after each line indicate which algorithm steps
    /// are being performed.
    ///
    /// [Renes-Costello-Batina 2015]: https://eprint.iacr.org/2015/1060
    pub fn add(&self, other: &ProjectivePoint) -> ProjectivePoint {
        debug_assert_eq!(
            CURVE_EQUATION_A,
            -FieldElement::from(3),
            "this implementation is only valid for C::EQUATION_A = -3"
        );

        let xx = self.x * other.x; // 1
        let yy = self.y * other.y; // 2
        let zz = self.z * other.z; // 3
        let xy_pairs = ((self.x + self.y) * (other.x + other.y)) - (xx + yy); // 4, 5, 6, 7, 8
        let yz_pairs = ((self.y + self.z) * (other.y + other.z)) - (yy + zz); // 9, 10, 11, 12, 13
        let xz_pairs = ((self.x + self.z) * (other.x + other.z)) - (xx + zz); // 14, 15, 16, 17, 18

        if cfg!(all(target_os = "zkvm", target_arch = "riscv32")) {
            let bzz_part = xz_pairs - zz.mul(CURVE_EQUATION_B); // 19, 20
            let bzz3_part = bzz_part.mul_single(3); // 21, 22

            let yy_m_bzz3 = yy - bzz3_part; // 23
            let yy_p_bzz3 = yy + bzz3_part; // 24

            let zz3 = zz.mul_single(3); // 26, 27
            let bxz_part = xz_pairs.mul(CURVE_EQUATION_B) - (zz3 + xx); // 25, 28, 29
            let bxz3_part = bxz_part.mul_single(3); // 30, 31
            let xx3_m_zz3 = xx.mul_single(3) - zz3; // 32, 33, 34

            return ProjectivePoint {
                x: (yy_p_bzz3.mul(xy_pairs)) - (yz_pairs.mul(bxz3_part)), // 35, 39, 40
                y: (yy_p_bzz3.mul(yy_m_bzz3)) + (xx3_m_zz3.mul(bxz3_part)), // 36, 37, 38
                z: (yy_m_bzz3.mul(yz_pairs)) + (xy_pairs.mul(xx3_m_zz3)), // 41, 42, 43
            };
        }

        let bzz_part = xz_pairs - (CURVE_EQUATION_B * zz); // 19, 20
        let bzz3_part = bzz_part.double() + bzz_part; // 21, 22
        let yy_m_bzz3 = yy - bzz3_part; // 23
        let yy_p_bzz3 = yy + bzz3_part; // 24

        let zz3 = zz.double() + zz; // 26, 27
        let bxz_part = (CURVE_EQUATION_B * xz_pairs) - (zz3 + xx); // 25, 28, 29
        let bxz3_part = bxz_part.double() + bxz_part; // 30, 31
        let xx3_m_zz3 = xx.double() + xx - zz3; // 32, 33, 34

        ProjectivePoint {
            x: (yy_p_bzz3 * xy_pairs) - (yz_pairs * bxz3_part), // 35, 39, 40
            y: (yy_p_bzz3 * yy_m_bzz3) + (xx3_m_zz3 * bxz3_part), // 36, 37, 38
            z: (yy_m_bzz3 * yz_pairs) + (xy_pairs * xx3_m_zz3), // 41, 42, 43
        }
    }

    /// Returns `self + other`.
    /// Implements complete mixed addition for curves with `a = -3`
    ///
    /// Implements the complete mixed addition formula from [Renes-Costello-Batina 2015]
    /// (Algorithm 5). The comments after each line indicate which algorithm
    /// steps are being performed.
    ///
    /// [Renes-Costello-Batina 2015]: https://eprint.iacr.org/2015/1060
    fn add_mixed(&self, other: &AffinePoint) -> ProjectivePoint {
        debug_assert_eq!(
            CURVE_EQUATION_A,
            -FieldElement::from(3),
            "this implementation is only valid for C::EQUATION_A = -3"
        );

        let xx = self.x * other.x; // 1
        let yy = self.y * other.y; // 2
        let xy_pairs = ((self.x + self.y) * (other.x + other.y)) - (xx + yy); // 3, 4, 5, 6, 7
        let yz_pairs = (other.y * self.z) + self.y; // 8, 9 (t4)
        let xz_pairs = (other.x * self.z) + self.x; // 10, 11 (y3)

        if cfg!(all(target_os = "zkvm", target_arch = "riscv32")) {
            let bz_part = xz_pairs - self.z.mul(CURVE_EQUATION_B); // 12, 13
            let bz3_part = bz_part.mul_single(3); // 14, 15
            let yy_m_bzz3 = yy - bz3_part; // 16
            let yy_p_bzz3 = yy + bz3_part; // 17

            let z3 = self.z.mul_single(3); // 19, 20
            let bxz_part = xz_pairs.mul(CURVE_EQUATION_B) - (z3 + xx); // 18, 21, 22
            let bxz3_part = bxz_part.mul_single(3); // 23, 24
            let xx3_m_zz3 = xx.mul_single(3) - z3; // 25, 26, 27

            let mut ret = ProjectivePoint {
                x: (yy_p_bzz3.mul(xy_pairs)) - (yz_pairs.mul(bxz3_part)), // 28, 32, 33
                y: (yy_p_bzz3.mul(yy_m_bzz3)) + (xx3_m_zz3.mul(bxz3_part)), // 29, 30, 31
                z: (yy_m_bzz3.mul(yz_pairs)) + (xy_pairs.mul(xx3_m_zz3)), // 34, 35, 36
            };
            ret.conditional_assign(self, other.is_identity());
            return ret;
        }

        let bz_part = xz_pairs - (CURVE_EQUATION_B * self.z); // 12, 13
        let bz3_part = bz_part.double() + bz_part; // 14, 15
        let yy_m_bzz3 = yy - bz3_part; // 16
        let yy_p_bzz3 = yy + bz3_part; // 17

        let z3 = self.z.double() + self.z; // 19, 20
        let bxz_part = (CURVE_EQUATION_B * xz_pairs) - (z3 + xx); // 18, 21, 22
        let bxz3_part = bxz_part.double() + bxz_part; // 23, 24
        let xx3_m_zz3 = xx.double() + xx - z3; // 25, 26, 27

        let mut ret = ProjectivePoint {
            x: (yy_p_bzz3 * xy_pairs) - (yz_pairs * bxz3_part), // 28, 32, 33
            y: (yy_p_bzz3 * yy_m_bzz3) + (xx3_m_zz3 * bxz3_part), // 29, 30, 31
            z: (yy_m_bzz3 * yz_pairs) + (xy_pairs * xx3_m_zz3), // 34, 35, 36
        };
        ret.conditional_assign(self, other.is_identity());
        ret
    }

    /// Returns `self - other`.
    pub fn sub(&self, other: &Self) -> Self {
        self.add(&other.neg())
    }

    /// Returns `self - other`.
    fn sub_mixed(&self, other: &AffinePoint) -> Self {
        self.add_mixed(&other.neg())
    }

    /// Returns `[k] self`.
    fn mul(&self, k: &Scalar) -> Self {
        // Into::into(*k) -> Uint for NIST P256 is U256
        let k = Into::<U256>::into(*k).to_le_byte_array();

        let mut pc = [Self::default(); 16];
        pc[0] = Self::IDENTITY;
        pc[1] = *self;

        for i in 2..16 {
            pc[i] = if i % 2 == 0 {
                Double::double(&pc[i / 2])
            } else {
                pc[i - 1].add(self)
            };
        }

        let mut q = Self::IDENTITY;
        let mut pos = U256::BITS - 4;

        loop {
            let slot = (k[pos >> 3] >> (pos & 7)) & 0xf;

            let mut t = ProjectivePoint::IDENTITY;

            for i in 1..16 {
                t.conditional_assign(
                    &pc[i],
                    Choice::from(((slot as usize ^ i).wrapping_sub(1) >> 8) as u8 & 1),
                );
            }

            q = q.add(&t);

            if pos == 0 {
                break;
            }

            q = Double::double(&Double::double(&Double::double(&Double::double(&q))));
            pos -= 4;
        }

        q
    }
}

impl CofactorGroup for ProjectivePoint {
    type Subgroup = Self;

    fn clear_cofactor(&self) -> Self::Subgroup {
        *self
    }

    fn into_subgroup(self) -> CtOption<Self> {
        CtOption::new(self, 1.into())
    }

    fn is_torsion_free(&self) -> Choice {
        1.into()
    }
}

impl ConditionallySelectable for ProjectivePoint {
    #[inline(always)]
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        Self {
            x: FieldElement::conditional_select(&a.x, &b.x, choice),
            y: FieldElement::conditional_select(&a.y, &b.y, choice),
            z: FieldElement::conditional_select(&a.z, &b.z, choice),
        }
    }
}

impl ConstantTimeEq for ProjectivePoint {
    fn ct_eq(&self, other: &Self) -> Choice {
        self.to_affine().ct_eq(&other.to_affine())
    }
}

impl Default for ProjectivePoint {
    fn default() -> Self {
        Self::IDENTITY
    }
}

impl DefaultIsZeroes for ProjectivePoint {}

impl Double for ProjectivePoint {
    /// Implements point doubling for curves with `a = -3`
    ///
    /// Implements the exception-free point doubling formula from [Renes-Costello-Batina 2015]
    /// (Algorithm 6). The comments after each line indicate which algorithm
    /// steps are being performed.
    ///
    /// [Renes-Costello-Batina 2015]: https://eprint.iacr.org/2015/1060
    fn double(&self) -> ProjectivePoint {
        debug_assert_eq!(
            CURVE_EQUATION_A,
            -FieldElement::from(3),
            "this implementation is only valid for C::EQUATION_A = -3"
        );

        let xx = self.x.square(); // 1
        let yy = self.y.square(); // 2
        let zz = self.z.square(); // 3
        let xy2 = (self.x * self.y).double(); // 4, 5
        let xz2 = (self.x * self.z).double(); // 6, 7

        if cfg!(all(target_os = "zkvm", target_arch = "riscv32")) {
            let bzz_part = zz.mul(CURVE_EQUATION_B) - xz2; // 8, 9
            let bzz3_part = bzz_part.mul_single(3); // 10, 11
            let yy_m_bzz3 = yy - bzz3_part; // 12
            let yy_p_bzz3 = yy + bzz3_part; // 13
            let y_frag = yy_p_bzz3.mul(yy_m_bzz3); // 14
            let x_frag = yy_m_bzz3.mul(xy2); // 15

            let zz3 = zz.mul_single(3); // 16, 17
            let bxz2_part = xz2.mul(CURVE_EQUATION_B) - (zz3 + xx); // 18, 19, 20
            let bxz6_part = bxz2_part.mul_single(3); // 21, 22
            let xx3_m_zz3 = xx.mul_single(3) - zz3; // 23, 24, 25

            let y = y_frag + (xx3_m_zz3.mul(bxz6_part)); // 26, 27
            let yz2 = (self.y.mul(self.z)).double();
            let x = x_frag - (bxz6_part.mul(yz2)); // 30, 31
            let z = (yz2.mul(yy)).mul_single(4); // 32, 33, 34

            return ProjectivePoint { x, y, z };
        }

        let bzz_part = (CURVE_EQUATION_B * zz) - xz2; // 8, 9
        let bzz3_part = bzz_part.double() + bzz_part; // 10, 11
        let yy_m_bzz3 = yy - bzz3_part; // 12
        let yy_p_bzz3 = yy + bzz3_part; // 13
        let y_frag = yy_p_bzz3 * yy_m_bzz3; // 14
        let x_frag = yy_m_bzz3 * xy2; // 15

        let zz3 = zz.double() + zz; // 16, 17
        let bxz2_part = (CURVE_EQUATION_B * xz2) - (zz3 + xx); // 18, 19, 20
        let bxz6_part = bxz2_part.double() + bxz2_part; // 21, 22
        let xx3_m_zz3 = xx.double() + xx - zz3; // 23, 24, 25

        let y = y_frag + (xx3_m_zz3 * bxz6_part); // 26, 27
        let yz2 = (self.y * self.z).double(); // 28, 29
        let x = x_frag - (bxz6_part * yz2); // 30, 31
        let z = (yz2 * yy).double().double(); // 32, 33, 34

        ProjectivePoint { x, y, z }
    }
}

impl Eq for ProjectivePoint {}

impl From<AffinePoint> for ProjectivePoint {
    fn from(p: AffinePoint) -> Self {
        let projective = ProjectivePoint {
            x: p.x,
            y: p.y,
            z: FieldElement::ONE,
        };
        Self::conditional_select(&projective, &Self::IDENTITY, p.is_identity())
    }
}

impl From<&AffinePoint> for ProjectivePoint {
    fn from(p: &AffinePoint) -> Self {
        Self::from(*p)
    }
}

impl From<PublicKey> for ProjectivePoint {
    fn from(public_key: PublicKey) -> ProjectivePoint {
        AffinePoint::from(public_key).into()
    }
}

impl From<&PublicKey> for ProjectivePoint {
    fn from(public_key: &PublicKey) -> ProjectivePoint {
        AffinePoint::from(public_key).into()
    }
}

impl FromEncodedPoint<NistP256> for ProjectivePoint {
    fn from_encoded_point(p: &EncodedPoint) -> CtOption<Self> {
        AffinePoint::from_encoded_point(p).map(Self::from)
    }
}

impl Group for ProjectivePoint {
    type Scalar = Scalar;

    fn random(mut rng: impl RngCore) -> Self {
        Self::GENERATOR * <Scalar as Field>::random(&mut rng)
    }

    fn identity() -> Self {
        Self::IDENTITY
    }

    fn generator() -> Self {
        Self::GENERATOR
    }

    fn is_identity(&self) -> Choice {
        self.ct_eq(&Self::IDENTITY)
    }

    #[must_use]
    fn double(&self) -> Self {
        Double::double(self)
    }
}

impl GroupEncoding for ProjectivePoint {
    type Repr = CompressedPoint;

    fn from_bytes(bytes: &Self::Repr) -> CtOption<Self> {
        <AffinePoint as GroupEncoding>::from_bytes(bytes).map(Into::into)
    }

    fn from_bytes_unchecked(bytes: &Self::Repr) -> CtOption<Self> {
        // No unchecked conversion possible for compressed points
        Self::from_bytes(bytes)
    }

    fn to_bytes(&self) -> Self::Repr {
        self.to_affine().to_bytes()
    }
}

impl group::Curve for ProjectivePoint {
    type AffineRepr = AffinePoint;

    fn to_affine(&self) -> AffinePoint {
        ProjectivePoint::to_affine(self)
    }

    // TODO(tarcieri): re-enable when we can add `Invert` bounds on `FieldElement`
    // #[cfg(feature = "alloc")]
    // #[inline]
    // fn batch_normalize(projective: &[Self], affine: &mut [Self::AffineRepr]) {
    //     assert_eq!(projective.len(), affine.len());
    //     let mut zs = vec![C::FieldElement::ONE; projective.len()];
    //     batch_normalize_generic(projective, zs.as_mut_slice(), affine);
    // }
}

impl<const N: usize> BatchNormalize<[ProjectivePoint; N]> for ProjectivePoint {
    type Output = [Self::AffineRepr; N];

    #[inline]
    fn batch_normalize(points: &[Self; N]) -> [Self::AffineRepr; N] {
        let mut zs = [FieldElement::ONE; N];
        let mut affine_points = [AffinePoint::IDENTITY; N];
        batch_normalize_generic(points, &mut zs, &mut affine_points);
        affine_points
    }
}

#[cfg(feature = "alloc")]
impl BatchNormalize<[ProjectivePoint]> for ProjectivePoint {
    type Output = Vec<Self::AffineRepr>;

    #[inline]
    fn batch_normalize(points: &[Self]) -> Vec<Self::AffineRepr> {
        let mut zs = vec![FieldElement::ONE; points.len()];
        let mut affine_points = vec![AffinePoint::IDENTITY; points.len()];
        batch_normalize_generic(points, zs.as_mut_slice(), &mut affine_points);
        affine_points
    }
}

// Generic implementation of batch normalization.
fn batch_normalize_generic<P, Z, O>(points: &P, zs: &mut Z, out: &mut O)
where
    FieldElement: BatchInvert<Z>,
    P: AsRef<[ProjectivePoint]> + ?Sized,
    Z: AsMut<[FieldElement]> + ?Sized,
    O: AsMut<[AffinePoint]> + ?Sized,
{
    let points = points.as_ref();
    let out = out.as_mut();

    for i in 0..points.len() {
        // Even a single zero value will fail inversion for the entire batch.
        // Put a dummy value (above `FieldElement::ONE`) so inversion succeeds
        // and treat that case specially later-on.
        zs.as_mut()[i].conditional_assign(&points[i].z, !points[i].z.ct_eq(&FieldElement::ZERO));
    }

    // This is safe to unwrap since we assured that all elements are non-zero
    let zs_inverses = <FieldElement as BatchInvert<Z>>::batch_invert(zs).unwrap();

    for i in 0..out.len() {
        // If the `z` coordinate is non-zero, we can use it to invert;
        // otherwise it defaults to the `IDENTITY` value.
        out[i] = AffinePoint::conditional_select(
            &points[i].to_affine_internal(zs_inverses.as_ref()[i]),
            &AffinePoint::IDENTITY,
            points[i].z.ct_eq(&FieldElement::ZERO),
        );
    }
}

impl LinearCombination for ProjectivePoint {}

impl MulByGenerator for ProjectivePoint {
    fn mul_by_generator(scalar: &Self::Scalar) -> Self {
        // TODO(tarcieri): precomputed basepoint tables
        Self::generator() * scalar
    }
}

impl PrimeGroup for ProjectivePoint
where
    Self: Double,
    FieldBytes<NistP256>: Copy,
    FieldBytesSize<NistP256>: ModulusSize,
    CompressedPoint: Copy,
    <UncompressedPointSize<NistP256> as ArrayLength<u8>>::ArrayType: Copy,
{
}

impl PrimeCurve for ProjectivePoint {
    type Affine = AffinePoint;
}

impl PartialEq for ProjectivePoint {
    fn eq(&self, other: &Self) -> bool {
        self.ct_eq(other).into()
    }
}

impl ToEncodedPoint<NistP256> for ProjectivePoint {
    fn to_encoded_point(&self, compress: bool) -> EncodedPoint {
        self.to_affine().to_encoded_point(compress)
    }
}

impl TryFrom<ProjectivePoint> for PublicKey {
    type Error = Error;

    fn try_from(point: ProjectivePoint) -> Result<PublicKey> {
        AffinePoint::from(point).try_into()
    }
}

impl TryFrom<&ProjectivePoint> for PublicKey {
    type Error = Error;

    fn try_from(point: &ProjectivePoint) -> Result<PublicKey> {
        AffinePoint::from(point).try_into()
    }
}

//
// Arithmetic trait impls
//

impl Add<ProjectivePoint> for ProjectivePoint {
    type Output = ProjectivePoint;

    fn add(self, other: ProjectivePoint) -> ProjectivePoint {
        ProjectivePoint::add(&self, &other)
    }
}

impl Add<&ProjectivePoint> for &ProjectivePoint {
    type Output = ProjectivePoint;

    fn add(self, other: &ProjectivePoint) -> ProjectivePoint {
        ProjectivePoint::add(self, other)
    }
}

impl Add<&ProjectivePoint> for ProjectivePoint {
    type Output = ProjectivePoint;

    fn add(self, other: &ProjectivePoint) -> ProjectivePoint {
        ProjectivePoint::add(&self, other)
    }
}

impl AddAssign<ProjectivePoint> for ProjectivePoint {
    fn add_assign(&mut self, rhs: ProjectivePoint) {
        *self = ProjectivePoint::add(self, &rhs);
    }
}

impl AddAssign<&ProjectivePoint> for ProjectivePoint {
    fn add_assign(&mut self, rhs: &ProjectivePoint) {
        *self = ProjectivePoint::add(self, rhs);
    }
}

impl Add<AffinePoint> for ProjectivePoint {
    type Output = ProjectivePoint;

    fn add(self, other: AffinePoint) -> ProjectivePoint {
        ProjectivePoint::add_mixed(&self, &other)
    }
}

impl Add<&AffinePoint> for &ProjectivePoint {
    type Output = ProjectivePoint;

    fn add(self, other: &AffinePoint) -> ProjectivePoint {
        ProjectivePoint::add_mixed(self, other)
    }
}

impl Add<&AffinePoint> for ProjectivePoint {
    type Output = ProjectivePoint;

    fn add(self, other: &AffinePoint) -> ProjectivePoint {
        ProjectivePoint::add_mixed(&self, other)
    }
}

impl AddAssign<AffinePoint> for ProjectivePoint {
    fn add_assign(&mut self, rhs: AffinePoint) {
        *self = ProjectivePoint::add_mixed(self, &rhs);
    }
}

impl AddAssign<&AffinePoint> for ProjectivePoint {
    fn add_assign(&mut self, rhs: &AffinePoint) {
        *self = ProjectivePoint::add_mixed(self, rhs);
    }
}

impl Sum for ProjectivePoint {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(ProjectivePoint::IDENTITY, |a, b| a + b)
    }
}

impl<'a> Sum<&'a ProjectivePoint> for ProjectivePoint {
    fn sum<I: Iterator<Item = &'a ProjectivePoint>>(iter: I) -> Self {
        iter.cloned().sum()
    }
}

impl Sub<ProjectivePoint> for ProjectivePoint {
    type Output = ProjectivePoint;

    fn sub(self, other: ProjectivePoint) -> ProjectivePoint {
        ProjectivePoint::sub(&self, &other)
    }
}

impl Sub<&ProjectivePoint> for &ProjectivePoint {
    type Output = ProjectivePoint;

    fn sub(self, other: &ProjectivePoint) -> ProjectivePoint {
        ProjectivePoint::sub(self, other)
    }
}

impl Sub<&ProjectivePoint> for ProjectivePoint {
    type Output = ProjectivePoint;

    fn sub(self, other: &ProjectivePoint) -> ProjectivePoint {
        ProjectivePoint::sub(&self, other)
    }
}

impl SubAssign<ProjectivePoint> for ProjectivePoint {
    fn sub_assign(&mut self, rhs: ProjectivePoint) {
        *self = ProjectivePoint::sub(self, &rhs);
    }
}

impl SubAssign<&ProjectivePoint> for ProjectivePoint {
    fn sub_assign(&mut self, rhs: &ProjectivePoint) {
        *self = ProjectivePoint::sub(self, rhs);
    }
}

impl Sub<AffinePoint> for ProjectivePoint {
    type Output = ProjectivePoint;

    fn sub(self, other: AffinePoint) -> ProjectivePoint {
        ProjectivePoint::sub_mixed(&self, &other)
    }
}

impl Sub<&AffinePoint> for &ProjectivePoint {
    type Output = ProjectivePoint;

    fn sub(self, other: &AffinePoint) -> ProjectivePoint {
        ProjectivePoint::sub_mixed(self, other)
    }
}

impl Sub<&AffinePoint> for ProjectivePoint {
    type Output = ProjectivePoint;

    fn sub(self, other: &AffinePoint) -> ProjectivePoint {
        ProjectivePoint::sub_mixed(&self, other)
    }
}

impl SubAssign<AffinePoint> for ProjectivePoint {
    fn sub_assign(&mut self, rhs: AffinePoint) {
        *self = ProjectivePoint::sub_mixed(self, &rhs);
    }
}

impl SubAssign<&AffinePoint> for ProjectivePoint {
    fn sub_assign(&mut self, rhs: &AffinePoint) {
        *self = ProjectivePoint::sub_mixed(self, rhs);
    }
}

impl Mul<Scalar> for ProjectivePoint {
    type Output = Self;

    fn mul(self, scalar: Scalar) -> Self {
        ProjectivePoint::mul(&self, scalar.borrow())
    }
}

impl Mul<&Scalar> for ProjectivePoint {
    type Output = ProjectivePoint;

    fn mul(self, scalar: &Scalar) -> ProjectivePoint {
        ProjectivePoint::mul(&self, &scalar)
    }
}

impl Mul<&Scalar> for &ProjectivePoint {
    type Output = ProjectivePoint;

    fn mul(self, scalar: &Scalar) -> ProjectivePoint {
        ProjectivePoint::mul(self, scalar)
    }
}

impl MulAssign<Scalar> for ProjectivePoint {
    fn mul_assign(&mut self, scalar: Scalar) {
        *self = ProjectivePoint::mul(self, scalar.borrow());
    }
}

impl MulAssign<&Scalar> for ProjectivePoint {
    fn mul_assign(&mut self, scalar: &Scalar) {
        *self = ProjectivePoint::mul(self, scalar);
    }
}

impl Neg for ProjectivePoint {
    type Output = ProjectivePoint;

    fn neg(self) -> ProjectivePoint {
        ProjectivePoint::neg(&self)
    }
}

impl<'a> Neg for &'a ProjectivePoint {
    type Output = ProjectivePoint;

    fn neg(self) -> ProjectivePoint {
        ProjectivePoint::neg(self)
    }
}
