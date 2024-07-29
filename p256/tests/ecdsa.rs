//! ECDSA tests.

#![cfg(feature = "arithmetic")]

use elliptic_curve::ops::Reduce;
use p256::{
    ecdsa::{SigningKey, VerifyingKey},
    NonZeroScalar, U256,
};
use proptest::prelude::*;

prop_compose! {
    fn signing_key()(bytes in any::<[u8; 32]>()) -> SigningKey {
        <NonZeroScalar as Reduce<U256>>::reduce_bytes(&bytes.into()).into()
    }
}

fn config() -> ProptestConfig {
    if cfg!(all(target_os = "zkvm", target_arch = "riscv32")) {
        ProptestConfig::with_cases(1)
    } else {
        ProptestConfig::default()
    }
}

proptest! {
    #![proptest_config(config())]

    #[test]
    fn recover_from_msg(sk in signing_key()) {
        let msg = b"example";
        let (signature, v) = sk.sign_recoverable(msg).unwrap();
        let recovered_vk = VerifyingKey::recover_from_msg(msg, &signature, v).unwrap();
        prop_assert_eq!(sk.verifying_key(), &recovered_vk);
    }
}
