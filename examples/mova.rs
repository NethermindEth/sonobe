use std::time::Instant;
use ark_ff::{PrimeField};
use folding_schemes::commitment::CommitmentScheme;
use folding_schemes::folding::mova::{Witness};
use folding_schemes::folding::mova::nifs::NIFS;
use folding_schemes::transcript::poseidon::{poseidon_canonical_config, PoseidonTranscript};
use folding_schemes::commitment::pedersen::{Pedersen};
use folding_schemes::utils::vec::{dense_matrix_to_sparse, SparseMatrix};
use ark_pallas::{Fr, Projective};
use ark_std::UniformRand;
use std::mem::size_of_val;
use folding_schemes::transcript::Transcript;
use ark_std::{log2};
use rand::Rng;
use folding_schemes::ccs::r1cs::R1CS;

fn main() {
    println!("starting");

    // define r1cs and parameters
    let random_num: usize = rand::thread_rng().gen_range(1..=2642245);
    let r1cs = get_test_r1cs();
    let mut rng = ark_std::test_rng();
    let (pedersen_params, _) = Pedersen::<Projective>::setup(&mut rng, r1cs.A.n_cols).unwrap();
    let poseidon_config = poseidon_canonical_config::<ark_pallas::Fr>();
    let mut transcript_p = PoseidonTranscript::<Projective>::new(&poseidon_config);

    // INSTANCE 1
    let z_1 = get_test_z(random_num);
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
    let result = NIFS::<Projective, Pedersen<Projective>,PoseidonTranscript<Projective>>::prove(&pedersen_params, &r1cs, &mut transcript_p, &committed_instance_1, &committed_instance_2, &witness_1, &witness_2);

    println!("Z number: {:?}", random_num);
    println!("Mova prove time {:?}", start.elapsed());
    println!("Mova bytes used {:?}", size_of_val(&result));

}


// would be best to move this to other file
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


