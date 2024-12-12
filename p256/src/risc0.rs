use risc0_bigint2::ec::{Curve, WeierstrassCurve, EC_256_WIDTH_WORDS};

/// The secp256r1 curve's prime field characteristic
pub(crate) const SECP256R1_PRIME: [u32; EC_256_WIDTH_WORDS] = [
    0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0x00000000, 0x00000000, 0x00000000, 0x00000001, 0xFFFFFFFF,
];

// /// The secp256r1 curve's order
// pub(crate) const SECP256R1_ORDER: [u32; EC_256_WIDTH_WORDS] = [
//     0xFC632551, 0xF3B9CAC2, 0xA7179E84, 0xBCE6FAAD, 0xFFFFFFFF, 0xFFFFFFFF, 0x00000000, 0xFFFFFFFF,
// ];

const SECP256R1_CURVE: &WeierstrassCurve<EC_256_WIDTH_WORDS> =
    &WeierstrassCurve::<EC_256_WIDTH_WORDS>::new(
        SECP256R1_PRIME,
        // Curve parameter a = -3 (represented mod p)
        [
            0xFFFFFFFC, 0xFFFFFFFF, 0xFFFFFFFF, 0x00000000, 0x00000000, 0x00000000, 0x00000001,
            0xFFFFFFFF,
        ],
        // Curve parameter b
        [
            0x27D2604B, 0x3BCE3C3E, 0xCC53B0F6, 0x651D06B0, 0x769886BC, 0xB3EBBD55, 0xAA3A93E7,
            0x5AC635D8,
        ],
    );


impl Curve<EC_256_WIDTH_WORDS> for crate::NistP256 {
    const CURVE: &'static WeierstrassCurve<EC_256_WIDTH_WORDS> = SECP256R1_CURVE;
}
