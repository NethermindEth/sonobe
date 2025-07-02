use ark_crypto_primitives::sponge::{
    merlin::Transcript as MerlinTranscript, Absorb, CryptographicSponge,
};
use ark_ec::{AffineRepr, CurveGroup};
use ark_ff::{BigInteger, PrimeField};

use super::{AbsorbNonNative, Transcript};

impl<F: PrimeField + Absorb> Transcript<F> for MerlinTranscript {
    fn absorb_point<C: CurveGroup<BaseField = F>>(&mut self, p: &C) {
        let (x, y) = p.into_affine().xy().unwrap_or_default();
        self.absorb(&x);
        self.absorb(&y);
    }
    fn absorb_nonnative<V: AbsorbNonNative>(&mut self, v: &V) {
        self.absorb(&v.to_native_sponge_field_elements_as_vec::<F>());
    }
    fn get_challenge(&mut self) -> F {
        let c = self.squeeze_field_elements(1);
        self.absorb(&c[0]);
        c[0]
    }
    fn get_challenge_nbits(&mut self, nbits: usize) -> Vec<bool> {
        let bits = self.squeeze_bits(nbits);
        self.absorb(&F::from(F::BigInt::from_bits_le(&bits)));
        bits
    }
    fn get_challenges(&mut self, n: usize) -> Vec<F> {
        let c = self.squeeze_field_elements(n);
        self.absorb(&c);
        c
    }
}
