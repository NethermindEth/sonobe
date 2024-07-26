/// Implements the scheme described in [Nova](https://eprint.iacr.org/2021/370.pdf) and
/// [CycleFold](https://eprint.iacr.org/2023/1192.pdf).
use ark_crypto_primitives::{
    crh::{poseidon::CRH, CRHScheme},
    sponge::{poseidon::PoseidonConfig, Absorb},
};
use ark_ec::{AffineRepr, CurveGroup, Group};
use ark_ff::{PrimeField, ToConstraintField};
use ark_poly::MultilinearExtension;

use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{fmt::Debug, log2};
use ark_std::{One, Zero};

use std::usize;

use crate::commitment::CommitmentScheme;

use crate::utils::vec::is_zero_vec;
use crate::Error;
use crate::transcript::{AbsorbNonNative, Transcript};

use crate::utils::mle::dense_vec_to_dense_mle;

pub mod homogenization;
pub mod nifs;
pub mod traits;
pub mod utils;

#[derive(Debug, Clone, Eq, PartialEq, CanonicalSerialize, CanonicalDeserialize)]
pub struct CommittedInstance<C: CurveGroup> {
    // Random evaluation point for the E
    pub rE: Vec<C::ScalarField>,
    // MLE of E
    pub mleE: C::ScalarField,
    pub u: C::ScalarField,
    pub cmW: C,
    pub x: Vec<C::ScalarField>,
}

impl<C: CurveGroup> CommittedInstance<C> {
    pub fn dummy(io_len: usize) -> Self {
        Self {
            rE: vec![C::ScalarField::zero(); io_len],
            mleE: C::ScalarField::zero(),
            u: C::ScalarField::zero(),
            cmW: C::zero(),
            x: vec![C::ScalarField::zero(); io_len],
        }
    }
}

impl<C: CurveGroup> Absorb for CommittedInstance<C>
    where
        C::ScalarField: Absorb,
{
    fn to_sponge_bytes(&self, dest: &mut Vec<u8>) {
        // This is never called
        unimplemented!()
    }

    fn to_sponge_field_elements<F: PrimeField>(&self, dest: &mut Vec<F>) {
        self.u.to_sponge_field_elements(dest);
        self.x.to_sponge_field_elements(dest);
        self.cmW
            .to_native_sponge_field_elements_as_vec()
            .to_sponge_field_elements(dest);
    }
}

impl<C: CurveGroup> AbsorbNonNative<C::BaseField> for CommittedInstance<C>
    where
        <C as ark_ec::CurveGroup>::BaseField: ark_ff::PrimeField + Absorb,
{
    fn to_native_sponge_field_elements(&self, dest: &mut Vec<C::BaseField>) {
        self.rE.to_native_sponge_field_elements(dest);
        [self.mleE].to_native_sponge_field_elements(dest);
        [self.u].to_native_sponge_field_elements(dest);
        self.x.to_native_sponge_field_elements(dest);
        let (cmW_x, cmW_y) = match self.cmW.into_affine().xy() {
            Some((&x, &y)) => (x, y),
            None => (C::BaseField::zero(), C::BaseField::zero()),
        };

        cmW_x.to_sponge_field_elements(dest);
        cmW_y.to_sponge_field_elements(dest);
    }
}

impl<C: CurveGroup> CommittedInstance<C>
    where
        <C as Group>::ScalarField: Absorb,
        <C as ark_ec::CurveGroup>::BaseField: ark_ff::PrimeField,
{
    /// hash implements the committed instance hash compatible with the gadget implemented in
    /// nova/circuits.rs::CommittedInstanceVar.hash.
    /// Returns `H(i, z_0, z_i, U_i)`, where `i` can be `i` but also `i+1`, and `U_i` is the
    /// `CommittedInstance`.
    pub fn hash<T: Transcript<C::ScalarField>>(
        &self,
        sponge: &T,
        pp_hash: C::ScalarField, // public params hash
        i: C::ScalarField,
        z_0: Vec<C::ScalarField>,
        z_i: Vec<C::ScalarField>,
    ) -> C::ScalarField {
        let mut sponge = sponge.clone();
        sponge.absorb(&pp_hash);
        sponge.absorb(&i);
        sponge.absorb(&z_0);
        sponge.absorb(&z_i);
        sponge.absorb(&self);
        sponge.squeeze_field_elements(1)[0]
    }
}



#[derive(Debug, Clone, Eq, PartialEq, CanonicalSerialize, CanonicalDeserialize)]
pub struct Witness<C: CurveGroup> {
    pub E: Vec<C::ScalarField>,
    pub W: Vec<C::ScalarField>,
    pub rW: C::ScalarField,
}

#[derive(Debug, Clone, Eq, PartialEq, CanonicalSerialize, CanonicalDeserialize)]
pub struct InstanceWitness<C: CurveGroup> {
    pub ci: CommittedInstance<C>,
    pub w: Witness<C>,
}

impl<C: CurveGroup> Witness<C>
where
    <C as Group>::ScalarField: Absorb,
{
    pub fn new(w: Vec<C::ScalarField>, e_len: usize) -> Self {
        // note: at the current version, we don't use the blinding factors and we set them to 0
        // always.
        Self {
            E: vec![C::ScalarField::zero(); e_len],
            W: w,
            rW: C::ScalarField::zero(),
        }
    }
    pub fn commit<CS: CommitmentScheme<C>>(
        &self,
        params: &CS::ProverParams,
        x: Vec<C::ScalarField>,
        rE: Vec<C::ScalarField>,
    ) -> Result<CommittedInstance<C>, Error> {
        let mut mleE = C::ScalarField::zero();
        if !is_zero_vec::<C::ScalarField>(&self.E) {
            let E = dense_vec_to_dense_mle(log2(self.E.len()) as usize, &self.E);
            mleE = E.evaluate(&rE).ok_or(Error::NotExpectedLength(
                rE.len(),
                log2(self.E.len()) as usize,
            ))?;
        }
        let cmW = CS::commit(params, &self.W, &self.rW)?;
        Ok(CommittedInstance {
            rE,
            mleE,
            u: C::ScalarField::one(),
            cmW,
            x,
        })
    }
}

