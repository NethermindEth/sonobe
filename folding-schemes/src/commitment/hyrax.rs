use ark_ec::{AffineRepr, CurveGroup};
use ark_ff::{Field, Zero};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::rand::RngCore;
use ark_std::vec::Vec;
use ark_std::{cfg_chunks, cfg_into_iter, UniformRand};
use rayon::prelude::*;
use std::fmt::Debug;

use crate::commitment::CommitmentScheme;
use crate::transcript::Transcript;
use crate::{Curve, Error};

#[derive(Debug, Clone, Eq, PartialEq, CanonicalSerialize, CanonicalDeserialize)]
pub struct Proof<C: Curve> {
    pub evaluation: Vec<C::ScalarField>,
}

#[derive(Clone, CanonicalSerialize, CanonicalDeserialize, Debug)]
pub struct HyraxGenerators<C: CurveGroup> {
    pub pedersen_generators: Vec<C::Affine>,
    pub row_len: usize,
    pub col_len: usize,
}

#[derive(Clone, CanonicalSerialize, CanonicalDeserialize, Debug)]
pub struct Hyrax<C: CurveGroup> {
    pub row_commitments: Vec<C>,
}

impl<C: CurveGroup> HyraxGenerators<C> {
    pub fn setup(len: usize) -> Self {
        let (col_len, row_len) = matrix_dimensions(len, 1);
        let pedersen_generators = (0..row_len).map(|_| C::Affine::generator()).collect();
        Self {
            pedersen_generators,
            row_len,
            col_len,
        }
    }
}

impl<C: CurveGroup> Hyrax<C> {
    pub fn commit(coeffs: &[C::ScalarField], gens: &HyraxGenerators<C>) -> Result<Self, Error> {
        assert_eq!(coeffs.len(), gens.row_len * gens.col_len);

        let row_commitments = cfg_chunks!(coeffs, gens.row_len)
            .map(|row| {
                let msm = C::msm_unchecked(&gens.pedersen_generators, row);
                Ok(msm)
            })
            .collect::<Result<Vec<_>, Error>>()?;
        Ok(Self { row_commitments })
    }

    pub fn open(
        coeffs: &[C::ScalarField],
        point: &[C::ScalarField],
        gens: &HyraxGenerators<C>,
    ) -> Result<Vec<C::ScalarField>, Error> {
        let l_vars = point[..gens.col_len.trailing_zeros() as usize].to_vec();
        let l_poly = eq_evals(&l_vars);
        let row_len = gens.row_len;
        let col_len = gens.col_len;

        let result = cfg_into_iter!(0..col_len)
            .map(|i| {
                let weight = l_poly[i];
                let offset = i * row_len;
                let mut contrib = vec![C::ScalarField::zero(); row_len];
                for j in 0..row_len {
                    contrib[j] = weight * coeffs[offset + j];
                }
                contrib
            })
            .reduce(
                || vec![C::ScalarField::zero(); row_len],
                |mut acc, contrib| {
                    for j in 0..row_len {
                        acc[j] += contrib[j];
                    }
                    acc
                },
            );

        Ok(result)
    }
}

fn eq_evals<F: Field>(point: &[F]) -> Vec<F> {
    let n = 1 << point.len();
    cfg_into_iter!(0..n)
        .map(|i| {
            let mut eval = F::one();
            for (j, &bit) in point.iter().enumerate() {
                let b = (i >> j) & 1 == 1;
                eval *= if b { bit } else { F::one() - bit };
            }
            eval
        })
        .collect()
}

fn matrix_dimensions(num_vars: usize, _ratio: usize) -> (usize, usize) {
    let row = (num_vars / 2).next_power_of_two();
    let col = (1 << num_vars) / row;
    (col, row)
}

impl<C: Curve, const H: bool> CommitmentScheme<C, H> for Hyrax<C> {
    type ProverParams = HyraxGenerators<C>;
    type VerifierParams = HyraxGenerators<C>;
    type Proof = Proof<C>;
    type ProverChallenge = ();
    type Challenge = ();

    fn is_hiding() -> bool {
        false
    }

    fn setup(
        mut rng: impl RngCore,
        len: usize,
    ) -> Result<(Self::ProverParams, Self::VerifierParams), Error> {
        let (col_len, row_len) = matrix_dimensions(len, 1);
        let pedersen_generators = (0..row_len).map(|_| C::Affine::rand(&mut rng)).collect();
        let params = HyraxGenerators {
            pedersen_generators,
            row_len,
            col_len,
        };
        Ok((params.clone(), params))
    }

    fn commit(
        params: &Self::ProverParams,
        v: &[C::ScalarField],
        _r: &C::ScalarField,
    ) -> Result<C, Error> {
        if v.len() != params.row_len * params.col_len {
            return Err(Error::PedersenParamsLen(
                params.row_len * params.col_len,
                v.len(),
            ));
        }

        let commitment = Hyrax::commit(v, params)?;
        Ok(commitment
            .row_commitments
            .iter()
            .copied()
            .reduce(|a, b| a + b)
            .unwrap_or_else(C::zero))
    }

    fn prove(
        params: &Self::ProverParams,
        _transcript: &mut impl Transcript<C::ScalarField>,
        _cm: &C,
        v: &[C::ScalarField],
        _r: &C::ScalarField,
        _rng: Option<&mut dyn RngCore>,
    ) -> Result<Self::Proof, Error> {
        let point = vec![C::ScalarField::zero(); (params.col_len as f64).log2() as usize];
        let evaluation = Hyrax::open(v, &point, params)?;
        Ok(Proof { evaluation })
    }

    fn prove_with_challenge(
        _params: &Self::ProverParams,
        _challenge: Self::ProverChallenge,
        _v: &[C::ScalarField],
        _r: &C::ScalarField,
        _rng: Option<&mut dyn RngCore>,
    ) -> Result<Self::Proof, Error> {
        unreachable!("Hyrax does not use prove_with_challenge")
    }

    fn verify(
        _params: &Self::VerifierParams,
        _transcript: &mut impl Transcript<C::ScalarField>,
        _cm: &C,
        _proof: &Self::Proof,
    ) -> Result<(), Error> {
        todo!("Verification not yet implemented")
    }

    fn verify_with_challenge(
        _params: &Self::VerifierParams,
        _challenge: Self::Challenge,
        _cm: &C,
        _proof: &Self::Proof,
    ) -> Result<(), Error> {
        unreachable!("Hyrax does not use verify_with_challenge")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::commitment::CommitmentScheme;
    use crate::transcript::poseidon::poseidon_canonical_config;
    use ark_crypto_primitives::sponge::poseidon::PoseidonSponge;
    use ark_crypto_primitives::sponge::CryptographicSponge;
    use ark_pallas::{Fr, Projective};
    use ark_std::{test_rng, UniformRand, Zero};

    #[test]
    fn test_hyrax_basic_commit_open() -> Result<(), Error> {
        let num_vars = 6;
        let (col_len, row_len) = matrix_dimensions(num_vars, 1);
        let total_len = col_len * row_len;

        let gens = HyraxGenerators::<Projective>::setup(num_vars);
        let mut rng = test_rng();
        let coeffs: Vec<Fr> = (0..total_len).map(|_| Fr::rand(&mut rng)).collect();
        let point: Vec<Fr> = (0..num_vars).map(|_| Fr::rand(&mut rng)).collect();

        let commitment = Hyrax::commit(&coeffs, &gens)?;
        assert_eq!(commitment.row_commitments.len(), col_len);

        let opening = Hyrax::open(&coeffs, &point, &gens)?;
        assert_eq!(opening.len(), row_len);
        Ok(())
    }

    #[test]
    fn test_hyrax_commitment_scheme_trait() -> Result<(), Error> {
        let num_vars = 6;
        let (params, _) =
            <Hyrax<Projective> as CommitmentScheme<Projective>>::setup(test_rng(), num_vars)?;

        let len = params.row_len * params.col_len;
        let mut rng = test_rng();
        let v: Vec<Fr> = (0..len).map(|_| Fr::rand(&mut rng)).collect();
        let r = Fr::zero();

        let cm = <Hyrax<Projective> as CommitmentScheme<Projective>>::commit(&params, &v, &r)?;

        let mut transcript = PoseidonSponge::new(&poseidon_canonical_config::<Fr>());
        let proof = <Hyrax<Projective> as CommitmentScheme<Projective>>::prove(
            &params,
            &mut transcript,
            &cm,
            &v,
            &r,
            None,
        )?;
        assert_eq!(proof.evaluation.len(), params.row_len);

        Ok(())
    }
}
