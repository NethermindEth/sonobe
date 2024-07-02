use ark_ff::PrimeField;
use ark_pallas::{Fr, Projective};
use ark_std::log2;
use ark_std::UniformRand;
use folding_schemes::ccs::r1cs::R1CS;
use folding_schemes::commitment::pedersen::Pedersen;
use folding_schemes::commitment::CommitmentScheme;
use folding_schemes::folding::mova::homogenization::{
    PointVsLineHomogenization, SumCheckHomogenization,
};
use folding_schemes::folding::mova::nifs::NIFS;
use folding_schemes::folding::mova::Witness;
use folding_schemes::transcript::poseidon::{poseidon_canonical_config, PoseidonTranscript};
use folding_schemes::transcript::Transcript;
use folding_schemes::utils::vec::{dense_matrix_to_sparse, SparseMatrix};
use num_bigint::{BigUint, RandBigInt};
use rand::Rng;
use std::mem::size_of_val;
use std::time::Instant;

use folding_schemes::Error;
use num_traits::{One, Zero};
use folding_schemes::folding::mova::traits::MovaR1CS;

fn main() {
    println!("starting");

    // define r1cs and parameters
    let r1cs: R1CS<Fr> = get_test_r1cs();
    let mut rng = ark_std::test_rng();
    let (pedersen_params, _) = Pedersen::<Projective>::setup(&mut rng, r1cs.A.n_cols).unwrap();
    let poseidon_config = poseidon_canonical_config::<ark_pallas::Fr>();
    let mut transcript_p = PoseidonTranscript::<Projective>::new(&poseidon_config);

    // let big_number: BigUint = One::one() ; // This creates a 250-bit number.

    // // INSTANCE 1
    let z_1: Vec<Fr> = get_test_z(3);
    let (w_1, x_1) = r1cs.split_z(&z_1);
    let witness_1 = Witness::<Projective>::new(w_1.clone(), r1cs.A.n_rows);

    // generate a random evaluation point for MLE
    let size_rE_1 = log2(witness_1.E.len());
    let rE_1: Vec<_> = (0..size_rE_1).map(|_| Fr::rand(&mut rng)).collect();

    let committed_instance_1 = witness_1
        .commit::<Pedersen<Projective>>(&pedersen_params, x_1, rE_1)
        .unwrap();

    // INSTANCE 2
    let z_2 = get_test_z(4);
    let (w_2, x_2) = r1cs.split_z(&z_2);
    let witness_2 = Witness::<Projective>::new(w_2.clone(), r1cs.A.n_rows);

    let size_rE_2 = log2(witness_2.E.len());
    let rE_2: Vec<_> = (0..size_rE_2).map(|_| Fr::rand(&mut rng)).collect();

    let mut committed_instance_2 = witness_2
        .commit::<Pedersen<Projective>>(&pedersen_params, x_2, rE_2)
        .unwrap();

    let start = Instant::now();
    // NIFS.P
    let result = NIFS::<
        Projective,
        Pedersen<Projective>,
        PoseidonTranscript<Projective>,
        PointVsLineHomogenization<Projective, PoseidonTranscript<Projective>>,
    >::prove(
        &pedersen_params,
        &r1cs,
        &mut transcript_p,
        &committed_instance_1,
        &committed_instance_2,
        &witness_1,
        &witness_2,
    ).unwrap();

    println!("Mova prove time {:?}", start.elapsed());
    println!("Mova bytes used {:?}", size_of_val(&result));

    //NIFS.V
    let poseidon_config = poseidon_canonical_config::<ark_pallas::Fr>();
    let mut transcript_p = PoseidonTranscript::<Projective>::new(&poseidon_config);
    let (proof, _, folded_w) = result;

    let folded_committed_instance = NIFS::<
        Projective,
        Pedersen<Projective>,
        PoseidonTranscript<Projective>,
        PointVsLineHomogenization<Projective, PoseidonTranscript<Projective>>,
    >::verify(&mut transcript_p, &committed_instance_1, &committed_instance_2, &proof).unwrap();
    let check = r1cs.check_relaxed_instance_relation(&folded_w, &folded_committed_instance);
    match check {
        Ok(_) => println!("The relation check was successful."),
        Err(e) => println!("The relation check failed: {:?}", e),
    }


}

pub fn get_test_r1cs<F: PrimeField>() -> R1CS<F> {
    // R1CS for: x^3 + x + 5 = y (example from article
    // https://www.vitalik.ca/general/2016/12/10/qap.html )
    let A = to_F_matrix::<F>(vec![
        vec![0, 1, 0, 0, 0, 0],
        vec![0, 0, 0, 1, 0, 0],
        vec![0, 1, 0, 0, 1, 0],
        vec![5, 0, 0, 0, 0, 1],
    ]);
    let B = to_F_matrix::<F>(vec![
        vec![0, 1, 0, 0, 0, 0],
        vec![0, 1, 0, 0, 0, 0],
        vec![1, 0, 0, 0, 0, 0],
        vec![1, 0, 0, 0, 0, 0],
    ]);
    let C = to_F_matrix::<F>(vec![
        vec![0, 0, 0, 1, 0, 0],
        vec![0, 0, 0, 0, 1, 0],
        vec![0, 0, 0, 0, 0, 1],
        vec![0, 0, 1, 0, 0, 0],
    ]);

    R1CS::<F> { l: 1, A, B, C }
}

pub fn get_test_z<F: PrimeField>(input: usize) -> Vec<F> {
    // z = (1, io, w)
    to_F_vec(vec![
        1,
        input,                             // io
        input * input * input + input + 5, // x^3 + x + 5
        input * input,                     // x^2
        input * input * input,             // x^2 * x
        input * input * input + input,     // x^3 + x
    ])
}

pub fn to_F_matrix<F: PrimeField>(M: Vec<Vec<usize>>) -> SparseMatrix<F> {
    dense_matrix_to_sparse(to_F_dense_matrix(M))
}
pub fn to_F_dense_matrix<F: PrimeField>(M: Vec<Vec<usize>>) -> Vec<Vec<F>> {
    M.iter()
        .map(|m| m.iter().map(|r| F::from(*r as u64)).collect())
        .collect()
}
pub fn to_F_vec<F: PrimeField>(z: Vec<usize>) -> Vec<F> {
    z.iter().map(|c| F::from(*c as u64)).collect()
}
