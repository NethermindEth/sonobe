use crate::bench_utils::{get_test_r1cs, get_test_z, write_to_csv};
use ark_pallas::{Fr, Projective};
use folding_schemes::commitment::pedersen::Pedersen;
use folding_schemes::commitment::CommitmentScheme;
use folding_schemes::transcript::poseidon::poseidon_canonical_config;
use std::time::{Duration, Instant};

use ark_crypto_primitives::sponge::poseidon::PoseidonSponge;
use ark_crypto_primitives::sponge::CryptographicSponge;
use ark_ff::PrimeField;
use ark_std::UniformRand;
use folding_schemes::arith::Arith;
use folding_schemes::arith::ccs::CCS;
use folding_schemes::arith::r1cs::R1CS;
use folding_schemes::folding::hypernova::nimfs::NIMFS;
use folding_schemes::utils::vec::SparseMatrix;

mod bench_utils;

fn hypernova_benchmarks(power: usize, prove_times: &mut Vec<Duration>) {
    // MOVA BENCH LINE
    // let r1cs: R1CS<Fr> = get_test_r1cs(power);
    let mut rng = ark_std::test_rng();
    // MOVA BENCH LINE
    // let ccs = CCS::<Fr>::from(r1cs);
    // NORMAL BENCH LINE
    let ccs = get_test_ccs::<Fr>();
    let (pedersen_params, _) = Pedersen::<Projective>::setup(&mut rng, ccs.n - ccs.l - 1).unwrap();
    // MOVA BENCH LINES
    // let z_1 = get_test_z(power, 384);
    // let z_2 = get_test_z(power, 20);
    // NORMAL BENCH LINE
    let z_1 = get_test_z2(3);
    let z_2 = get_test_z2(4);


    let (running_instance, w1) = ccs
        .to_lcccs::<_, _, Pedersen<Projective>, false>(&mut rng, &pedersen_params, &z_1)
        .unwrap();
    // MOVA BENCH LINE
    // running_instance.u = Fr::rand(&mut rng);

    let (new_instance, w2) = ccs
        .to_cccs::<_, _, Pedersen<Projective>, false>(&mut rng, &pedersen_params, &z_2)
        .unwrap();


    let poseidon_config = poseidon_canonical_config::<Fr>();

    let mut transcript_p: PoseidonSponge<Fr> = PoseidonSponge::<Fr>::new(&poseidon_config);
    transcript_p.absorb(&Fr::from_le_bytes_mod_order(b"init init"));

    let start = Instant::now();

    let (proof, folded_lcccs, folded_witness, _) = NIMFS::<Projective, PoseidonSponge<Fr>>::prove(
        &mut transcript_p,
        &ccs,
        &[running_instance.clone()],
        &[new_instance.clone()],
        &[w1],
        &[w2],
    )
        .unwrap();

    let prove_time = start.elapsed();
    prove_times.push(prove_time);
    println!("Hypernova prove time {:?}", prove_time);

    // MOVA BENCH LINES -- THE REST SHOULD BE COMMENTED OUT
    let mut transcript_v: PoseidonSponge<Fr> = PoseidonSponge::<Fr>::new(&poseidon_config);
    transcript_v.absorb(&Fr::from_le_bytes_mod_order(b"init init"));
    //
    // // Run the verifier side of the multifolding
    let folded_lcccs_v = NIMFS::<Projective, PoseidonSponge<Fr>>::verify(
        &mut transcript_v,
        &ccs,
        &[running_instance.clone()],
        &[new_instance.clone()],
        proof,
    )
        .unwrap();
    assert_eq!(folded_lcccs, folded_lcccs_v);

    // Check that the folded LCCCS instance is a valid instance with respect to the folded witness
    ccs.check_relation(&folded_witness, &folded_lcccs).unwrap();
}

fn main() {
    // let pows: Vec<usize> = (10..24).collect();
    let pows: Vec<usize> = vec![10];
    let iter = 1;
    let mut prove_times: Vec<Duration> = Vec::with_capacity(pows.len() * iter);
    for i in 0..iter {
        println!("starting {:}", i);

        println!("{:?}", pows);

        for pow in &pows {
            println!("{}", pow);
            hypernova_benchmarks(*pow, &mut prove_times);
        }

        println!("Powers {:?}", pows);
        println!("Prove times {:?}", prove_times);
    }
    if let Err(e) = write_to_csv(&pows, &prove_times, "hypernova_prove_times.csv".to_string()) {
        eprintln!("Failed to write to CSV: {}", e);
    } else {
        println!("CSV file has been successfully written.");
    }
}

pub fn get_test_ccs<F: PrimeField>() -> CCS<F> {
    get_test_r1cs2::<F>().into()
}

pub fn get_test_r1cs2<F: PrimeField>() -> R1CS<F> {
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

pub fn get_test_z2<F: PrimeField>(input: usize) -> Vec<F> {
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

pub fn to_F_vec<F: PrimeField>(z: Vec<usize>) -> Vec<F> {
    z.iter().map(|c| F::from(*c as u64)).collect()
}

pub fn to_F_matrix<F: PrimeField>(M: Vec<Vec<usize>>) -> SparseMatrix<F> {
    dense_matrix_to_sparse(to_F_dense_matrix(M))
}
pub fn to_F_dense_matrix<F: PrimeField>(M: Vec<Vec<usize>>) -> Vec<Vec<F>> {
    M.iter()
        .map(|m| m.iter().map(|r| F::from(*r as u64)).collect())
        .collect()
}

pub fn dense_matrix_to_sparse<F: PrimeField>(m: Vec<Vec<F>>) -> SparseMatrix<F> {
    let mut r = SparseMatrix::<F> {
        n_rows: m.len(),
        n_cols: m[0].len(),
        coeffs: Vec::new(),
    };
    for m_row in m.iter() {
        let mut row: Vec<(F, usize)> = Vec::new();
        for (col_i, value) in m_row.iter().enumerate() {
            if !value.is_zero() {
                row.push((*value, col_i));
            }
        }
        r.coeffs.push(row);
    }
    r
}