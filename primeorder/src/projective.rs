//! Projective curve points.

#![allow(clippy::needless_range_loop, clippy::op_ref)]

use crate::{point_arithmetic::PointArithmetic, AffinePoint, Field, PrimeCurveParams};
use core::{
    borrow::Borrow,
    iter::Sum,
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};
use elliptic_curve::{
    bigint::{ArrayEncoding, Integer},
    generic_array::ArrayLength,
    group::{
        self,
        cofactor::CofactorGroup,
        prime::{PrimeCurve, PrimeGroup},
        Group, GroupEncoding,
    },
    ops::{LinearCombination, MulByGenerator},
    point::Double,
    rand_core::RngCore,
    sec1::{
        CompressedPoint, EncodedPoint, FromEncodedPoint, ModulusSize, ToEncodedPoint,
        UncompressedPointSize,
    },
    subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption},
    zeroize::DefaultIsZeroes,
    Error, FieldBytes, FieldBytesSize, PublicKey, Result, Scalar,
};

/// Point on a Weierstrass curve in projective coordinates.
#[derive(Clone, Copy, Debug)]
pub struct ProjectivePoint<C: PrimeCurveParams> {
    pub(crate) x: C::FieldElement,
    pub(crate) y: C::FieldElement,
    pub(crate) z: C::FieldElement,
}

impl<C> ProjectivePoint<C>
where
    C: PrimeCurveParams,
{
    /// Additive identity of the group a.k.a. the point at infinity.
    pub const IDENTITY: Self = Self {
        x: C::FieldElement::ZERO,
        y: C::FieldElement::ONE,
        z: C::FieldElement::ZERO,
    };

    /// Base point of the curve.
    pub const GENERATOR: Self = Self {
        x: C::GENERATOR.0,
        y: C::GENERATOR.1,
        z: C::FieldElement::ONE,
    };

    /// Returns the affine representation of this point, or `None` if it is the identity.
    pub fn to_affine(&self) -> AffinePoint<C> {
        self.z
            .invert()
            .map(|zinv| AffinePoint {
                x: self.x * &zinv,
                y: self.y * &zinv,
                infinity: 0,
            })
            .unwrap_or(AffinePoint::IDENTITY)
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
    pub fn add(&self, other: &Self) -> Self {
        C::PointArithmetic::add(self, other)
    }

    /// Returns `self + other`.
    fn add_mixed(&self, other: &AffinePoint<C>) -> Self {
        C::PointArithmetic::add_mixed(self, other)
    }

    /// Returns `self - other`.
    pub fn sub(&self, other: &Self) -> Self {
        self.add(&other.neg())
    }

    /// Returns `self - other`.
    fn sub_mixed(&self, other: &AffinePoint<C>) -> Self {
        self.add_mixed(&other.neg())
    }

    /// Returns `[k] self`.
    fn mul(&self, k: &Scalar<C>) -> Self
    where
        Self: Double,
    {
        #[cfg(all(target_os = "zkvm", target_arch = "riscv32"))]
       {
            let scalar = scalar_to_words::<C>(k);
            let affine = projective_to_affine(self);

            let mut result = ec::AffinePoint::new_unchecked([0u32; 8], [0u32; 8]);
            affine.mul(&scalar, &mut result);
            return affine_to_projective(&result);
        }
        #[cfg(not(all(target_os = "zkvm", target_arch = "riscv32")))]
        {
            let k = Into::<C::Uint>::into(*k).to_le_byte_array();

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
            let mut pos = C::Uint::BITS - 4;

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
}

#[cfg(all(target_os = "zkvm", target_arch = "riscv32"))]
use risc0_bigint2::ec;

#[cfg(all(target_os = "zkvm", target_arch = "riscv32"))]
use elliptic_curve::{generic_array::GenericArray, PrimeField};

#[cfg(all(target_os = "zkvm", target_arch = "riscv32"))]
// TODO remove inline
#[inline(never)]
fn projective_to_affine<C>(p: &ProjectivePoint<C>) -> ec::AffinePoint<8, C>
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

#[cfg(all(target_os = "zkvm", target_arch = "riscv32"))]
// TODO remove inline
#[inline(never)]
fn affine_to_projective<C>(affine: &ec::AffinePoint<8, C>) -> ProjectivePoint<C>
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

#[cfg(all(target_os = "zkvm", target_arch = "riscv32"))]
fn scalar_to_words<C>(s: &Scalar<C>) -> [u32; 8]
where
    C: PrimeCurveParams,
{
    let mut bytes: [u8; 32] = s.to_repr().as_slice().try_into().unwrap();
    // U256 is big endian, need to flip to little endian.
    bytes.reverse();
    bytemuck::cast::<_, [u32; 8]>(bytes)
}

impl<C> CofactorGroup for ProjectivePoint<C>
where
    Self: Double,
    C: PrimeCurveParams,
    FieldBytes<C>: Copy,
    FieldBytesSize<C>: ModulusSize,
    CompressedPoint<C>: Copy,
    <UncompressedPointSize<C> as ArrayLength<u8>>::ArrayType: Copy,
{
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

impl<C> ConditionallySelectable for ProjectivePoint<C>
where
    C: PrimeCurveParams,
{
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        Self {
            x: C::FieldElement::conditional_select(&a.x, &b.x, choice),
            y: C::FieldElement::conditional_select(&a.y, &b.y, choice),
            z: C::FieldElement::conditional_select(&a.z, &b.z, choice),
        }
    }
}

impl<C> ConstantTimeEq for ProjectivePoint<C>
where
    C: PrimeCurveParams,
{
    fn ct_eq(&self, other: &Self) -> Choice {
        self.to_affine().ct_eq(&other.to_affine())
    }
}

impl<C> Default for ProjectivePoint<C>
where
    C: PrimeCurveParams,
{
    fn default() -> Self {
        Self::IDENTITY
    }
}

impl<C> DefaultIsZeroes for ProjectivePoint<C> where C: PrimeCurveParams {}

impl<C: PrimeCurveParams> Double for ProjectivePoint<C> {
    fn double(&self) -> Self {
        C::PointArithmetic::double(self)
    }
}

impl<C> Eq for ProjectivePoint<C> where C: PrimeCurveParams {}

impl<C> From<AffinePoint<C>> for ProjectivePoint<C>
where
    C: PrimeCurveParams,
{
    fn from(p: AffinePoint<C>) -> Self {
        let projective = ProjectivePoint {
            x: p.x,
            y: p.y,
            z: C::FieldElement::ONE,
        };
        Self::conditional_select(&projective, &Self::IDENTITY, p.is_identity())
    }
}

impl<C> From<&AffinePoint<C>> for ProjectivePoint<C>
where
    C: PrimeCurveParams,
{
    fn from(p: &AffinePoint<C>) -> Self {
        Self::from(*p)
    }
}

impl<C> From<PublicKey<C>> for ProjectivePoint<C>
where
    C: PrimeCurveParams,
{
    fn from(public_key: PublicKey<C>) -> ProjectivePoint<C> {
        AffinePoint::from(public_key).into()
    }
}

impl<C> From<&PublicKey<C>> for ProjectivePoint<C>
where
    C: PrimeCurveParams,
{
    fn from(public_key: &PublicKey<C>) -> ProjectivePoint<C> {
        AffinePoint::<C>::from(public_key).into()
    }
}

impl<C> FromEncodedPoint<C> for ProjectivePoint<C>
where
    C: PrimeCurveParams,
    FieldBytes<C>: Copy,
    FieldBytesSize<C>: ModulusSize,
    CompressedPoint<C>: Copy,
{
    fn from_encoded_point(p: &EncodedPoint<C>) -> CtOption<Self> {
        AffinePoint::<C>::from_encoded_point(p).map(Self::from)
    }
}

impl<C> Group for ProjectivePoint<C>
where
    Self: Double,
    C: PrimeCurveParams,
{
    type Scalar = Scalar<C>;

    fn random(mut rng: impl RngCore) -> Self {
        Self::GENERATOR * <Scalar<C> as Field>::random(&mut rng)
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

impl<C> GroupEncoding for ProjectivePoint<C>
where
    C: PrimeCurveParams,
    FieldBytes<C>: Copy,
    FieldBytesSize<C>: ModulusSize,
    CompressedPoint<C>: Copy,
    <UncompressedPointSize<C> as ArrayLength<u8>>::ArrayType: Copy,
{
    type Repr = CompressedPoint<C>;

    fn from_bytes(bytes: &Self::Repr) -> CtOption<Self> {
        <AffinePoint<C> as GroupEncoding>::from_bytes(bytes).map(Into::into)
    }

    fn from_bytes_unchecked(bytes: &Self::Repr) -> CtOption<Self> {
        // No unchecked conversion possible for compressed points
        Self::from_bytes(bytes)
    }

    fn to_bytes(&self) -> Self::Repr {
        self.to_affine().to_bytes()
    }
}

impl<C> group::Curve for ProjectivePoint<C>
where
    Self: Double,
    C: PrimeCurveParams,
{
    type AffineRepr = AffinePoint<C>;

    fn to_affine(&self) -> AffinePoint<C> {
        ProjectivePoint::to_affine(self)
    }
}

impl<C> LinearCombination for ProjectivePoint<C>
where
    Self: Double,
    C: PrimeCurveParams,
{
}

impl<C> MulByGenerator for ProjectivePoint<C>
where
    Self: Double,
    C: PrimeCurveParams,
{
    fn mul_by_generator(scalar: &Self::Scalar) -> Self {
        // TODO(tarcieri): precomputed basepoint tables
        Self::generator() * scalar
    }
}

impl<C> PrimeGroup for ProjectivePoint<C>
where
    Self: Double,
    C: PrimeCurveParams,
    FieldBytes<C>: Copy,
    FieldBytesSize<C>: ModulusSize,
    CompressedPoint<C>: Copy,
    <UncompressedPointSize<C> as ArrayLength<u8>>::ArrayType: Copy,
{
}

impl<C> PrimeCurve for ProjectivePoint<C>
where
    Self: Double,
    C: PrimeCurveParams,
    FieldBytes<C>: Copy,
    FieldBytesSize<C>: ModulusSize,
    CompressedPoint<C>: Copy,
    <UncompressedPointSize<C> as ArrayLength<u8>>::ArrayType: Copy,
{
    type Affine = AffinePoint<C>;
}

impl<C> PartialEq for ProjectivePoint<C>
where
    C: PrimeCurveParams,
{
    fn eq(&self, other: &Self) -> bool {
        self.ct_eq(other).into()
    }
}

impl<C> ToEncodedPoint<C> for ProjectivePoint<C>
where
    C: PrimeCurveParams,
    FieldBytesSize<C>: ModulusSize,
    CompressedPoint<C>: Copy,
    <UncompressedPointSize<C> as ArrayLength<u8>>::ArrayType: Copy,
{
    fn to_encoded_point(&self, compress: bool) -> EncodedPoint<C> {
        self.to_affine().to_encoded_point(compress)
    }
}

impl<C> TryFrom<ProjectivePoint<C>> for PublicKey<C>
where
    C: PrimeCurveParams,
{
    type Error = Error;

    fn try_from(point: ProjectivePoint<C>) -> Result<PublicKey<C>> {
        AffinePoint::<C>::from(point).try_into()
    }
}

impl<C> TryFrom<&ProjectivePoint<C>> for PublicKey<C>
where
    C: PrimeCurveParams,
{
    type Error = Error;

    fn try_from(point: &ProjectivePoint<C>) -> Result<PublicKey<C>> {
        AffinePoint::<C>::from(point).try_into()
    }
}

//
// Arithmetic trait impls
//

impl<C> Add<ProjectivePoint<C>> for ProjectivePoint<C>
where
    C: PrimeCurveParams,
{
    type Output = ProjectivePoint<C>;

    fn add(self, other: ProjectivePoint<C>) -> ProjectivePoint<C> {
        ProjectivePoint::add(&self, &other)
    }
}

impl<C> Add<&ProjectivePoint<C>> for &ProjectivePoint<C>
where
    C: PrimeCurveParams,
{
    type Output = ProjectivePoint<C>;

    fn add(self, other: &ProjectivePoint<C>) -> ProjectivePoint<C> {
        ProjectivePoint::add(self, other)
    }
}

impl<C> Add<&ProjectivePoint<C>> for ProjectivePoint<C>
where
    C: PrimeCurveParams,
{
    type Output = ProjectivePoint<C>;

    fn add(self, other: &ProjectivePoint<C>) -> ProjectivePoint<C> {
        ProjectivePoint::add(&self, other)
    }
}

impl<C> AddAssign<ProjectivePoint<C>> for ProjectivePoint<C>
where
    C: PrimeCurveParams,
{
    fn add_assign(&mut self, rhs: ProjectivePoint<C>) {
        *self = ProjectivePoint::add(self, &rhs);
    }
}

impl<C> AddAssign<&ProjectivePoint<C>> for ProjectivePoint<C>
where
    C: PrimeCurveParams,
{
    fn add_assign(&mut self, rhs: &ProjectivePoint<C>) {
        *self = ProjectivePoint::add(self, rhs);
    }
}

impl<C> Add<AffinePoint<C>> for ProjectivePoint<C>
where
    C: PrimeCurveParams,
{
    type Output = ProjectivePoint<C>;

    fn add(self, other: AffinePoint<C>) -> ProjectivePoint<C> {
        ProjectivePoint::add_mixed(&self, &other)
    }
}

impl<C> Add<&AffinePoint<C>> for &ProjectivePoint<C>
where
    C: PrimeCurveParams,
{
    type Output = ProjectivePoint<C>;

    fn add(self, other: &AffinePoint<C>) -> ProjectivePoint<C> {
        ProjectivePoint::add_mixed(self, other)
    }
}

impl<C> Add<&AffinePoint<C>> for ProjectivePoint<C>
where
    C: PrimeCurveParams,
{
    type Output = ProjectivePoint<C>;

    fn add(self, other: &AffinePoint<C>) -> ProjectivePoint<C> {
        ProjectivePoint::add_mixed(&self, other)
    }
}

impl<C> AddAssign<AffinePoint<C>> for ProjectivePoint<C>
where
    C: PrimeCurveParams,
{
    fn add_assign(&mut self, rhs: AffinePoint<C>) {
        *self = ProjectivePoint::add_mixed(self, &rhs);
    }
}

impl<C> AddAssign<&AffinePoint<C>> for ProjectivePoint<C>
where
    C: PrimeCurveParams,
{
    fn add_assign(&mut self, rhs: &AffinePoint<C>) {
        *self = ProjectivePoint::add_mixed(self, rhs);
    }
}

impl<C> Sum for ProjectivePoint<C>
where
    C: PrimeCurveParams,
{
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(ProjectivePoint::IDENTITY, |a, b| a + b)
    }
}

impl<'a, C> Sum<&'a ProjectivePoint<C>> for ProjectivePoint<C>
where
    C: PrimeCurveParams,
{
    fn sum<I: Iterator<Item = &'a ProjectivePoint<C>>>(iter: I) -> Self {
        iter.cloned().sum()
    }
}

impl<C> Sub<ProjectivePoint<C>> for ProjectivePoint<C>
where
    C: PrimeCurveParams,
{
    type Output = ProjectivePoint<C>;

    fn sub(self, other: ProjectivePoint<C>) -> ProjectivePoint<C> {
        ProjectivePoint::sub(&self, &other)
    }
}

impl<C> Sub<&ProjectivePoint<C>> for &ProjectivePoint<C>
where
    C: PrimeCurveParams,
{
    type Output = ProjectivePoint<C>;

    fn sub(self, other: &ProjectivePoint<C>) -> ProjectivePoint<C> {
        ProjectivePoint::sub(self, other)
    }
}

impl<C> Sub<&ProjectivePoint<C>> for ProjectivePoint<C>
where
    C: PrimeCurveParams,
{
    type Output = ProjectivePoint<C>;

    fn sub(self, other: &ProjectivePoint<C>) -> ProjectivePoint<C> {
        ProjectivePoint::sub(&self, other)
    }
}

impl<C> SubAssign<ProjectivePoint<C>> for ProjectivePoint<C>
where
    C: PrimeCurveParams,
{
    fn sub_assign(&mut self, rhs: ProjectivePoint<C>) {
        *self = ProjectivePoint::sub(self, &rhs);
    }
}

impl<C> SubAssign<&ProjectivePoint<C>> for ProjectivePoint<C>
where
    C: PrimeCurveParams,
{
    fn sub_assign(&mut self, rhs: &ProjectivePoint<C>) {
        *self = ProjectivePoint::sub(self, rhs);
    }
}

impl<C> Sub<AffinePoint<C>> for ProjectivePoint<C>
where
    C: PrimeCurveParams,
{
    type Output = ProjectivePoint<C>;

    fn sub(self, other: AffinePoint<C>) -> ProjectivePoint<C> {
        ProjectivePoint::sub_mixed(&self, &other)
    }
}

impl<C> Sub<&AffinePoint<C>> for &ProjectivePoint<C>
where
    C: PrimeCurveParams,
{
    type Output = ProjectivePoint<C>;

    fn sub(self, other: &AffinePoint<C>) -> ProjectivePoint<C> {
        ProjectivePoint::sub_mixed(self, other)
    }
}

impl<C> Sub<&AffinePoint<C>> for ProjectivePoint<C>
where
    C: PrimeCurveParams,
{
    type Output = ProjectivePoint<C>;

    fn sub(self, other: &AffinePoint<C>) -> ProjectivePoint<C> {
        ProjectivePoint::sub_mixed(&self, other)
    }
}

impl<C> SubAssign<AffinePoint<C>> for ProjectivePoint<C>
where
    C: PrimeCurveParams,
{
    fn sub_assign(&mut self, rhs: AffinePoint<C>) {
        *self = ProjectivePoint::sub_mixed(self, &rhs);
    }
}

impl<C> SubAssign<&AffinePoint<C>> for ProjectivePoint<C>
where
    C: PrimeCurveParams,
{
    fn sub_assign(&mut self, rhs: &AffinePoint<C>) {
        *self = ProjectivePoint::sub_mixed(self, rhs);
    }
}

impl<C, S> Mul<S> for ProjectivePoint<C>
where
    Self: Double,
    C: PrimeCurveParams,
    S: Borrow<Scalar<C>>,
{
    type Output = Self;

    fn mul(self, scalar: S) -> Self {
        ProjectivePoint::mul(&self, scalar.borrow())
    }
}

impl<C> Mul<&Scalar<C>> for &ProjectivePoint<C>
where
    C: PrimeCurveParams,
    ProjectivePoint<C>: Double,
{
    type Output = ProjectivePoint<C>;

    fn mul(self, scalar: &Scalar<C>) -> ProjectivePoint<C> {
        ProjectivePoint::mul(self, scalar)
    }
}

impl<C, S> MulAssign<S> for ProjectivePoint<C>
where
    Self: Double,
    C: PrimeCurveParams,
    S: Borrow<Scalar<C>>,
{
    fn mul_assign(&mut self, scalar: S) {
        *self = ProjectivePoint::mul(self, scalar.borrow());
    }
}

impl<C> Neg for ProjectivePoint<C>
where
    C: PrimeCurveParams,
{
    type Output = ProjectivePoint<C>;

    fn neg(self) -> ProjectivePoint<C> {
        ProjectivePoint::neg(&self)
    }
}

impl<'a, C> Neg for &'a ProjectivePoint<C>
where
    C: PrimeCurveParams,
{
    type Output = ProjectivePoint<C>;

    fn neg(self) -> ProjectivePoint<C> {
        ProjectivePoint::neg(self)
    }
}
