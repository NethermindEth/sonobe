use ark_ec::AffineRepr;
use ark_ff::{Field, Zero};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::iterable::Iterable;
use ark_std::rand::RngCore;
use ark_std::vec::Vec;
use rayon::prelude::*;
use std::fmt::Debug;
use std::marker::PhantomData;

use crate::commitment::pedersen::{Params, Pedersen};
use crate::commitment::{CommitmentScheme, NethermindCommitmentScheme};
use crate::{Curve, Error};

/// Taken from jolt but we assume ratio is 1 since we are dealing with square matrices
fn matrix_dimensions(num_vars: usize) -> (usize, usize) {
    let left_num_vars = num_vars / 2;
    let right_num_vars = num_vars - left_num_vars;

    let col_size = 2_usize.pow(left_num_vars as u32);
    let row_size = 2_usize.pow(right_num_vars as u32);

    (col_size , row_size)
}

#[derive(Clone, CanonicalSerialize, CanonicalDeserialize, Debug)]
pub struct HyraxGenerators<C: Curve> {
    pub pedersen_generators: Params<C>,
}

#[derive(Clone, CanonicalSerialize, CanonicalDeserialize, Debug, PartialEq, Eq)]
pub struct Hyrax<C: Curve> {
    _c: PhantomData<C>,
}

impl<C: Curve> HyraxGenerators<C> {
    pub fn setup(rng: impl RngCore, len: usize) -> Self {
        let (_col_len, row_len) = matrix_dimensions(len);
        let gens = Pedersen::<C, false>::setup2(rng, row_len).unwrap();
        HyraxGenerators {
            pedersen_generators: gens,
        }
    }
}

impl<C: Curve> Hyrax<C> {
    pub fn commit(elems: &[C::ScalarField], gens: &HyraxGenerators<C>) -> Result<Vec<C>, Error> {
        let n = elems.len();
        let ell = {
            if n.is_power_of_two() {
                (1usize.leading_zeros() - n.leading_zeros()) as usize
            } else {
                (0usize.leading_zeros() - n.leading_zeros()) as usize
            }
        };

        let (L_size, R_size) = matrix_dimensions(ell);
        assert_eq!(L_size * R_size, n);

        let row_commitments: Vec<C> = elems
            .par_chunks(R_size)
            .map(|row| {
                Pedersen::<C, false>::commit(
                    &gens.pedersen_generators,
                    row,
                    &C::ScalarField::zero(),
                )
                .unwrap()
            })
            .collect();
        Ok(row_commitments)
    }

    pub fn batch_commit(
        batch: &[&[C::ScalarField]],
        gens: &HyraxGenerators<C>,
    ) -> Result<Vec<Vec<C>>, Error> {
        let n = batch[0].len();
        batch.iter().for_each(|poly| assert_eq!(poly.len(), n));
        let ell = {
            if n.is_power_of_two() {
                (1usize.leading_zeros() - n.leading_zeros()) as usize
            } else {
                (0usize.leading_zeros() - n.leading_zeros()) as usize
            }
        };
        let (L_size, R_size) = matrix_dimensions(ell);
        assert_eq!(L_size * R_size, n);

        let rows = batch.par_iter().flat_map(|poly| poly.par_chunks(R_size));
        let row_commitments: Vec<C> = rows
            .map(|row| {
                Pedersen::<C, false>::commit(
                    &gens.pedersen_generators,
                    row,
                    &C::ScalarField::zero(),
                )
                .unwrap()
            })
            .collect();

        Ok(row_commitments
            .par_chunks(L_size)
            .map(|chunk| chunk.to_vec())
            .collect())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_pallas::{Fr, Projective};
    use ark_std::{test_rng, UniformRand};

    #[test]
    fn test_matrix_dimensions() {
        assert_eq!(matrix_dimensions(0), (1, 1));
        assert_eq!(matrix_dimensions(1), (1, 2));
        assert_eq!(matrix_dimensions(2), (2, 2));
        assert_eq!(matrix_dimensions(4), (4, 4));
        assert_eq!(matrix_dimensions(6), (8, 8));

        let num_vars = 30; // Very large number of variables
        let (cols, rows) = matrix_dimensions(num_vars);

        let expected_left = num_vars / 2;
        let expected_right = num_vars - expected_left;

        assert_eq!(cols, 1 << expected_left);
        assert_eq!(rows, 1 << expected_right);
    }

    #[test]
    fn test_setup() {
        let mut rng = test_rng();
        let gens = HyraxGenerators::<Projective>::setup(&mut rng, 4);
        assert!(!gens.pedersen_generators.generators.is_empty());
    }

    #[test]
    fn test_commit_success() {
        let mut rng = test_rng();
        let len = 4; // Must be a power of two number of bits
        let elems: Vec<Fr> = (0..(1 << len)).map(|_| Fr::rand(&mut rng)).collect();

        let gens = HyraxGenerators::<Projective>::setup(&mut rng, len);
        let commitment = Hyrax::<Projective>::commit(&elems, &gens);
        assert!(commitment.is_ok());

        let hyrax = commitment.unwrap();
        let (l_size, _) = matrix_dimensions(len);
        assert_eq!(hyrax.len(), l_size);
    }

    #[test]
    fn test_same_commit() {
        let mut rng = test_rng();
        let len = 4; // Must be a power of two number of bits
        let elems: Vec<Fr> = (0..(1 << len)).map(|_| Fr::rand(&mut rng)).collect();

        let gens = HyraxGenerators::<Projective>::setup(&mut rng, len);
        let commitment = Hyrax::<Projective>::commit(&elems, &gens);
        assert!(commitment.is_ok());

        let hyrax = commitment.unwrap();

        let commitment2 = Hyrax::<Projective>::commit(&elems, &gens);
        assert!(commitment2.is_ok());

        let hyrax2 = commitment2.unwrap();
        assert_eq!(hyrax, hyrax2);
    }

    #[test]
    fn test_batch_commit() {
        let mut rng = test_rng();
        let len = 3;
        let poly_len = 1 << len;
        let batch_size = 5;

        let batch: Vec<Vec<Fr>> = (0..batch_size)
            .map(|_| (0..poly_len).map(|_| Fr::rand(&mut rng)).collect())
            .collect();

        let batch_refs: Vec<&[Fr]> = batch.iter().map(|v| v.as_slice()).collect();

        let gens = HyraxGenerators::<Projective>::setup(&mut rng, len);
        let result = Hyrax::<Projective>::batch_commit(&batch_refs, &gens);
        assert!(result.is_ok());

        let commitments = result.unwrap();
        assert_eq!(commitments.len(), batch_size);
        let (l_size, _) = matrix_dimensions(len);
        for c in commitments {
            assert_eq!(c.len(), l_size);
        }
    }
}
