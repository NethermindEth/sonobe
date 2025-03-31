use ark_ec::{AffineRepr, CurveGroup};
use ark_ff::{Field, Zero};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::rand::RngCore;
use ark_std::vec::Vec;
// use ark_std::{cfg_chunks, cfg_into_iter, UniformRand};
use rayon::prelude::*;
use std::fmt::Debug;
use std::io::Read;
use rand_chacha::ChaCha20Rng;


use crate::commitment::{CommitmentScheme, NethermindCommitmentScheme};
use crate::transcript::Transcript;
use crate::{Curve, Error};
use crate::commitment::pedersen::{Params, Pedersen};
use ark_std::log2;
use rand::SeedableRng;
use sha3::digest::{ExtendableOutput, Update};
use sha3::Shake256;


/// Taken from jolt but we assume ratio is 1 since we are dealing with square matrices
fn matrix_dimensions(num_vars: usize) -> (usize, usize) {
    let left_num_vars = num_vars / 2;
    let right_num_vars = num_vars - left_num_vars;

    let col_size = 2_i32.pow(left_num_vars as u32);
    let row_size = 2_i32.pow(right_num_vars as u32);

    (col_size as usize, row_size as usize)
}

#[derive(Clone, CanonicalSerialize, CanonicalDeserialize, Debug)]
pub struct HyraxGenerators<C: Curve> {
    pub pedersen_generators: Params<C>,
}

#[derive(Clone, CanonicalSerialize, CanonicalDeserialize, Debug)]
pub struct Hyrax<C: Curve> {
    pub row_commitments: Vec<C>,
}

// impl<C: Curve> HyraxGenerators<C> {
//     pub fn setup(mut rng: impl RngCore, len: usize) -> Self {
//         let (_col_len, row_len) = matrix_dimensions(len);
//         let gens = Pedersen::<C, false>::setup2(rng, row_len).unwrap();
//         HyraxGenerators{pedersen_generators: gens}
//
//     }
// }

impl<G: Curve> HyraxGenerators<G>{
    pub fn new(num_vars: usize) -> Self {
        let (left, right) = matrix_dimensions(num_vars);
        let gens = Self::new2(right, b"Jolt v1 Hyrax generators");
        let h = G::zero();
        let params = Params{h, generators:  CurveGroup::normalize_batch(&gens[..right])};
        HyraxGenerators{pedersen_generators: params}
    }

    pub fn new2(len: usize, label: &[u8]) -> Vec<G> {
        let mut shake = Shake256::default();
        shake.update(label);
        let mut buf = vec![];
        G::generator().serialize_compressed(&mut buf).unwrap();
        shake.update(&buf);

        let mut reader = shake.finalize_xof();
        let mut seed = [0u8; 32];
        reader.read_exact(&mut seed).unwrap();
        let mut rng = ChaCha20Rng::from_seed(seed);

        let mut generators: Vec<G> = Vec::new();
        for _ in 0..len {
            generators.push(G::rand(&mut rng));
        }
        println!("len gen {}", generators.len());

        generators
    }
}


impl<C: Curve> Hyrax<C> {
    pub fn commit(elems: &[C::ScalarField], gens: &HyraxGenerators<C>) -> Result<Self, Error> {
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

        // let gens = CurveGroup::normalize_batch(&gens.pedersen_generators.generators[..R_size]);
        let row_commitments: Vec<C> = elems
            .par_chunks(R_size)
            .map(|row| Pedersen::<C, false>::commit(&gens.pedersen_generators, row, &C::ScalarField::zero()).unwrap())
            .collect();
        Ok(Self { row_commitments })
    }

//     pub fn open(
//         coeffs: &[C::ScalarField],
//         point: &[C::ScalarField],
//         gens: &HyraxGenerators<C>,
//     ) -> Result<Vec<C::ScalarField>, Error> {
//         let l_vars = point[..gens.col_len.trailing_zeros() as usize].to_vec();
//         let l_poly = eq_evals(&l_vars);
//         let row_len = gens.row_len;
//         let col_len = gens.col_len;
//
//         let result = cfg_into_iter!(0..col_len)
//             .map(|i| {
//                 let weight = l_poly[i];
//                 let offset = i * row_len;
//                 let mut contrib = vec![C::ScalarField::zero(); row_len];
//                 for j in 0..row_len {
//                     contrib[j] = weight * coeffs[offset + j];
//                 }
//                 contrib
//             })
//             .reduce(
//                 || vec![C::ScalarField::zero(); row_len],
//                 |mut acc, contrib| {
//                     for j in 0..row_len {
//                         acc[j] += contrib[j];
//                     }
//                     acc
//                 },
//             );
//
//         Ok(result)
//     }
// }
//
// fn eq_evals<F: Field>(point: &[F]) -> Vec<F> {
//     let n = 1 << point.len();
//     cfg_into_iter!(0..n)
//         .map(|i| {
//             let mut eval = F::one();
//             for (j, &bit) in point.iter().enumerate() {
//                 let b = (i >> j) & 1 == 1;
//                 eval *= if b { bit } else { F::one() - bit };
//             }
//             eval
//         })
//         .collect()
}

// #[cfg(test)]
// mod tests {
//     use super::*;
//     use crate::commitment::CommitmentScheme;
//     use crate::transcript::poseidon::poseidon_canonical_config;
//     use ark_crypto_primitives::sponge::poseidon::PoseidonSponge;
//     use ark_crypto_primitives::sponge::CryptographicSponge;
//     use ark_pallas::{Fr, Projective};
//     use ark_std::{test_rng, UniformRand, Zero};
//
//     #[test]
//     fn test_hyrax_basic_commit_open() -> Result<(), Error> {
//         let num_vars = 6;
//         let (col_len, row_len) = matrix_dimensions(num_vars, 1);
//         let total_len = col_len * row_len;
//
//         let gens = HyraxGenerators::<Projective>::setup(num_vars);
//         let mut rng = test_rng();
//         let coeffs: Vec<Fr> = (0..total_len).map(|_| Fr::rand(&mut rng)).collect();
//         let point: Vec<Fr> = (0..num_vars).map(|_| Fr::rand(&mut rng)).collect();
//
//         let commitment = Hyrax::commit(&coeffs, &gens)?;
//         assert_eq!(commitment.row_commitments.len(), col_len);
//
//         let opening = Hyrax::open(&coeffs, &point, &gens)?;
//         assert_eq!(opening.len(), row_len);
//         Ok(())
//     }
//
//     #[test]
//     fn test_hyrax_commitment_scheme_trait() -> Result<(), Error> {
//         let num_vars = 6;
//         let (params, _) =
//             <Hyrax<Projective> as CommitmentScheme<Projective>>::setup(test_rng(), num_vars)?;
//
//         let len = params.row_len * params.col_len;
//         let mut rng = test_rng();
//         let v: Vec<Fr> = (0..len).map(|_| Fr::rand(&mut rng)).collect();
//         let r = Fr::zero();
//
//         let cm = <Hyrax<Projective> as CommitmentScheme<Projective>>::commit(&params, &v, &r)?;
//
//         let mut transcript = PoseidonSponge::new(&poseidon_canonical_config::<Fr>());
//         let proof = <Hyrax<Projective> as CommitmentScheme<Projective>>::prove(
//             &params,
//             &mut transcript,
//             &cm,
//             &v,
//             &r,
//             None,
//         )?;
//         assert_eq!(proof.evaluation.len(), params.row_len);
//
//         Ok(())
//     }
// }

#[cfg(test)]
mod tests {
    use ark_bn254::{Fr, G1Projective as G1};
    use ark_poly_commit::hyrax::HyraxCommitment;
    use ark_std::log2;
    use crate::commitment::hyrax::{Hyrax, HyraxGenerators};

    #[test]
    fn compare_with_jolt(){
        let mut rng = ark_std::test_rng();
        let num_vars = 16;
        let matrix: Vec<Fr> = (0..16).map(|i| Fr::from(i as u64)).collect();
        println!("matrix {:?}", matrix);

        let params = HyraxGenerators::<G1>::new(log2(num_vars) as usize);
        println!("params {:?}", params.pedersen_generators.generators);

        let commitment = Hyrax::commit(&matrix, &params);
        println!("commitment {:?}", commitment);


    }
}
