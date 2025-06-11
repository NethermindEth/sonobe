use crate::commitment::pedersen::{Params, Pedersen};
use crate::commitment::{CommitmentScheme, NethermindCommitmentScheme};
use crate::{Curve, Error};
use ark_ff::Zero;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::iterable::Iterable;
use ark_std::log2;
use ark_std::rand::RngCore;
use ark_std::vec::Vec;
use rayon::prelude::*;
use std::fmt::Debug;
use std::marker::PhantomData;

/// Taken from Jolt but we assume ratio is 1 since we are dealing with square matrices
fn matrix_dimensions(num_elems: usize) -> (usize, usize) {
    let num_vars = if num_elems == 1 {
        1
    } else {
        log2(num_elems) as usize
    };

    let mut row_size = 2_usize.pow((num_vars / 2) as u32);
    row_size = row_size.next_power_of_two();

    let right_num_vars: usize = std::cmp::min(log2(row_size) as usize, num_vars - 1);
    row_size = 2_usize.pow(right_num_vars as u32);
    let left_num_vars = num_vars - right_num_vars;
    let col_size = 2_usize.pow(left_num_vars as u32);

    (col_size, row_size)
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
        let gens = Pedersen::<C, false>::setup_prover(rng, row_len).unwrap();
        HyraxGenerators {
            pedersen_generators: gens,
        }
    }
}

impl<C: Curve> Hyrax<C> {
    pub fn commit(elems: &[C::ScalarField], gens: &HyraxGenerators<C>) -> Result<Vec<C>, Error> {
        let n = elems.len();

        let (L_size, R_size) = matrix_dimensions(n);
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

        let (L_size, R_size) = matrix_dimensions(n);
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

    pub fn commit_sparse_matrix(
        indices_values: &[(usize, C::ScalarField)],
        gens: &HyraxGenerators<C>,
    ) -> Result<Vec<C>, Error> {
        let max_elems =
            gens.pedersen_generators.generators.len() * gens.pedersen_generators.generators.len();
        let (L_size, R_size) = matrix_dimensions(max_elems);
        // For each row i in [0..L_size], gather all (pos, val) where row_start <= pos < row_end
        // and do a "sparse" Pedersen commit using the row‐local indices (pos - row_start).
        let row_commitments: Vec<C> = (0..L_size)
            .into_par_iter()
            .map(|row_idx| {
                let row_start = row_idx * R_size;
                let row_end = row_start + R_size;

                // gather the pairs that fall in row i
                let row_sparse: Vec<(usize, C::ScalarField)> = indices_values
                    .iter()
                    .filter_map(|(pos, val)| {
                        if *pos >= row_start && *pos < row_end {
                            // shift the index for the row commit
                            Some((pos - row_start, *val))
                        } else {
                            None
                        }
                    })
                    .collect();

                Pedersen::<C, false>::commit_sparse(
                    &gens.pedersen_generators,
                    &row_sparse,
                    &C::ScalarField::zero(),
                )
            })
            .collect::<Result<Vec<C>, _>>()?;

        Ok(row_commitments)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_pallas::{Fr, Projective};
    use ark_std::{test_rng, UniformRand};
    use rand::Rng;
    #[test]
    fn test_matrix_dimensions() {
        // Values taken from the original Jolt implementation
        assert_eq!(matrix_dimensions(1), (2, 1));
        assert_eq!(matrix_dimensions(4), (2, 2));
        assert_eq!(matrix_dimensions(8), (4, 2));
        assert_eq!(matrix_dimensions(16), (4, 4));
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
        let dim = 4; // Must be a power of two number of bits
        let elems: Vec<Fr> = (0..(1 << dim)).map(|_| Fr::rand(&mut rng)).collect();

        let gens = HyraxGenerators::<Projective>::setup(&mut rng, elems.len());
        let commitment = Hyrax::<Projective>::commit(&elems, &gens);
        assert!(commitment.is_ok());

        let hyrax = commitment.unwrap();
        let (l_size, _) = matrix_dimensions(elems.len());
        assert_eq!(hyrax.len(), l_size);
    }

    #[test]
    fn test_same_commit() {
        let mut rng = test_rng();
        let len = 4; // Must be a power of two
        let elems: Vec<Fr> = (0..(1 << len)).map(|_| Fr::rand(&mut rng)).collect();

        let gens = HyraxGenerators::<Projective>::setup(&mut rng, elems.len());
        let commitment = Hyrax::<Projective>::commit(&elems, &gens);
        assert!(commitment.is_ok());

        let hyrax = commitment.unwrap();

        let commitment2 = Hyrax::<Projective>::commit(&elems, &gens);
        assert!(commitment2.is_ok());

        let hyrax2 = commitment2.unwrap();
        assert_eq!(hyrax, hyrax2);
    }

    #[test]
    fn test_sparse_dense() {
        // Create a sparse matrix
        let mut rng = test_rng();
        let dim = 8;
        let max_elems = 1 << dim;
        let mut elems = vec![Fr::zero(); max_elems];
        let random_idx: Vec<usize> = (0..dim).map(|_| rng.gen_range(0..dim * dim)).collect();
        for idx in random_idx {
            elems[idx] = Fr::rand(&mut rng);
        }

        // Compute its dense commitment
        let gens = HyraxGenerators::<Projective>::setup(&mut rng, elems.len());
        let commitment = Hyrax::<Projective>::commit(&elems, &gens);
        assert!(commitment.is_ok());
        let hyrax = commitment.unwrap();

        // Compute its sparse commitment
        let sparse_repr = elems
            .into_iter()
            .enumerate()
            .filter_map(|(idx, elem)| {
                if !elem.is_zero() {
                    Some((idx, elem))
                } else {
                    None
                }
            })
            .collect::<Vec<_>>();
        let commitment2 = Hyrax::<Projective>::commit_sparse_matrix(&sparse_repr, &gens);
        assert!(commitment2.is_ok());
        let hyrax2 = commitment2.unwrap();

        // Make sure they are equal
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

        let gens = HyraxGenerators::<Projective>::setup(&mut rng, poly_len);
        let result = Hyrax::<Projective>::batch_commit(&batch_refs, &gens);
        assert!(result.is_ok());

        let commitments = result.unwrap();
        assert_eq!(commitments.len(), batch_size);
        let (l_size, _) = matrix_dimensions(poly_len);
        for c in commitments {
            assert_eq!(c.len(), l_size);
        }
    }
}
