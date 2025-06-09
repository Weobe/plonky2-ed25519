use num::{BigUint, Integer};
use plonky2::field::types::{Field, PrimeField};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha512};

use crate::curve::curve_types::{AffinePoint, Curve};
use crate::curve::ed25519::Ed25519;
use crate::curve::ed25519::mul_naive;
use crate::field::ed25519_base::Ed25519Base;
use crate::field::ed25519_scalar::Ed25519Scalar;

pub const SAMPLE_MSG1: &str = "test message";
pub const SAMPLE_MSG2: &str = "plonky2";
pub const SAMPLE_PK1: [u8; 32] = [
    59, 106, 39, 188, 206, 182, 164, 45, 98, 163, 168, 208, 42, 111, 13, 115, 101, 50, 21, 119, 29,
    226, 67, 166, 58, 192, 72, 161, 139, 89, 218, 41,
];
pub const SAMPLE_SIG1: [u8; 64] = [
    104, 196, 204, 44, 176, 120, 225, 128, 47, 67, 245, 210, 247, 65, 201, 66, 34, 159, 217, 32,
    175, 224, 14, 12, 31, 231, 83, 160, 214, 122, 250, 68, 250, 203, 33, 143, 184, 13, 247, 140,
    185, 25, 122, 25, 253, 195, 83, 102, 240, 255, 30, 21, 108, 249, 77, 184, 36, 72, 9, 198, 49,
    12, 68, 8,
];
pub const SAMPLE_SIG2: [u8; 64] = [
    130, 82, 60, 170, 184, 218, 199, 182, 66, 19, 182, 14, 141, 214, 229, 180, 43, 19, 227, 183,
    130, 204, 69, 112, 171, 113, 6, 111, 218, 227, 249, 85, 57, 216, 145, 63, 71, 192, 201, 10, 54,
    234, 203, 8, 63, 240, 226, 101, 84, 167, 36, 246, 153, 35, 31, 52, 244, 82, 239, 137, 18, 62,
    134, 7,
];

#[derive(Copy, Clone, Debug, Deserialize, Eq, Hash, PartialEq, Serialize)]
pub struct EDDSASignature<C: Curve> {
    pub r: AffinePoint<C>,
    pub s: C::ScalarField,
}

///  √‑1  mod p  (same value as `POWER_OF_TWO_GENERATOR`)
const SQRT_M1: Ed25519Base = Ed25519Base::POWER_OF_TWO_GENERATOR;

/// Decompress a 32‑byte Ed25519 point into your own `AffinePoint`.
pub fn point_decompress(bytes: &[u8]) -> Option<AffinePoint<Ed25519>> {
    if bytes.len() != 32 {
        return None;
    }

    // --- 1. split into sign bit + y  -----------------------------
    let mut buf = [0u8; 32];
    buf.copy_from_slice(bytes);
    let sign = (buf[31] >> 7) & 1;
    buf[31] &= 0x7F; // clear sign bit ⇒ y bytes
    let y_big = BigUint::from_bytes_le(&buf);

    if &y_big >= &Ed25519Base::order() {
        return None;
    }
    let y = Ed25519Base::from_noncanonical_biguint(y_big);

    // y^2
    let y2 = y * y;

    // --- 2. compute x² = (y² – 1)/(d y² + 1) mod p --------------
    let num = y2 - Ed25519Base::ONE;
    let denom = Ed25519::D * y2 + Ed25519Base::ONE;
    let x2 = num * denom.inverse(); // denom ≠ 0 because |y| < p

    // --- 3. take square root  -----------------------------------
    //  p ≡ 5 (mod 8) ⇒  x = x²^((p+3)/8)
    let e = (Ed25519Base::order() + BigUint::from(3u8)) >> 3;
    let mut x = x2.exp_biguint(&e);

    // If x² != x2, multiply by √‑1
    if x * x != x2 {
        x *= SQRT_M1;
    }
    if x * x != x2 {
        return None;
    } // not a square ⇒ invalid

    // --- 4. correct sign  ---------------------------------------
    let x_is_odd = x.to_canonical_biguint().bit(0) as u8; // bool → 0/1
    if x_is_odd != sign {
        x = -x;
    }

    Some(AffinePoint::nonzero(x, y))
}

pub fn verify_message(msg: &[u8], sigv: &[u8], pkv: &[u8]) -> bool {
    let mut data = Vec::new();
    data.extend_from_slice(&sigv[..32]);
    data.extend_from_slice(pkv);
    data.extend_from_slice(msg);
    let data_u8 = data.as_slice();

    let mut hasher = Sha512::new();
    hasher.update(data_u8);
    let hash = hasher.finalize();
    let h_big_int = BigUint::from_bytes_le(hash.as_slice());
    let h_mod_25519 = h_big_int.mod_floor(&Ed25519Scalar::order());
    let h = Ed25519Scalar::from_noncanonical_biguint(h_mod_25519);

    let pk = point_decompress(pkv).expect("invalid public key");
    assert!(pk.is_valid());
    let r = point_decompress(&sigv[..32]).expect("invalid R");
    assert!(r.is_valid());
    let s = Ed25519Scalar::from_noncanonical_biguint(BigUint::from_bytes_le(&sigv[32..]));

    let g = Ed25519::GENERATOR_PROJECTIVE;
    let sb = mul_naive(s, g);
    let ha = mul_naive(h, pk.to_projective());
    let rhs = r + ha.to_affine();

    sb.to_affine() == rhs
}

#[cfg(test)]
mod tests {
    use crate::curve::eddsa::{
        SAMPLE_MSG1, SAMPLE_MSG2, SAMPLE_PK1, SAMPLE_SIG1, SAMPLE_SIG2, verify_message,
    };

    #[test]
    fn test_eddsa_native() {
        let result = verify_message(
            SAMPLE_MSG1.as_bytes(),
            SAMPLE_SIG1.as_slice(),
            SAMPLE_PK1.as_slice(),
        );
        assert!(result);
        let result = verify_message(
            SAMPLE_MSG2.as_bytes(),
            SAMPLE_SIG2.as_slice(),
            SAMPLE_PK1.as_slice(),
        );
        assert!(result);
    }
}

#[cfg(test)]
mod dalek_roundtrip {
    use super::verify_message; // your verifier
    use ed25519_dalek::{Keypair, Signer}; // Signer::sign(msg)
    use rand_core::OsRng; // same rand_core 0.6 as dalek

    /// Convert a dalek signature (R||S) into raw 64‑byte array
    fn sig_bytes(sig: &ed25519_dalek::Signature) -> [u8; 64] {
        sig.to_bytes() // already R‖S in dalek‑1.x
    }

    #[test]
    fn roundtrip_with_dalek() {
        // 1. random key‑pair
        let kp: Keypair = Keypair::generate(&mut OsRng);
        let pk_bytes: [u8; 32] = kp.public.to_bytes();

        // 2. two demo messages
        let msg1 = b"dalek test message";
        let msg2 = b"another example";

        // 3. sign with dalek (Sha‑512 is baked in)
        let sig1 = kp.sign(msg1);
        let sig2 = kp.sign(msg2);

        // 4. verify with *your* code
        assert!(verify_message(msg1, &sig_bytes(&sig1), &pk_bytes));
        assert!(verify_message(msg2, &sig_bytes(&sig2), &pk_bytes));
    }
}
