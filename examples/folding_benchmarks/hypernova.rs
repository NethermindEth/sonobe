use crate::bench_utils::{get_test_r1cs, get_test_z, write_to_csv};
use ark_ff::{BigInteger, Field, PrimeField};
use ark_pallas::{Fr, Projective};
use ark_std::{log2, UniformRand};
use folding_schemes::commitment::pedersen::Pedersen;
use folding_schemes::commitment::CommitmentScheme;
use folding_schemes::transcript::poseidon::poseidon_canonical_config;
use folding_schemes::transcript::Transcript;
use folding_schemes::utils::sum_check::SumCheck;
use rand::Rng;
use std::mem::size_of_val;
use std::time::{Duration, Instant};

use ark_crypto_primitives::sponge::poseidon::PoseidonSponge;
use ark_crypto_primitives::sponge::CryptographicSponge;
use folding_schemes::arith::ccs::CCS;
use folding_schemes::arith::r1cs::R1CS;
use folding_schemes::folding::hypernova::nimfs::NIMFS;
use folding_schemes::utils::vec::{dense_matrix_to_sparse, SparseMatrix};
use std::error::Error;

mod bench_utils;

fn hypernova_benchmarks(power: usize, prove_times: &mut Vec<Duration>) {
    let size = 1 << power;

    let r1cs: R1CS<Fr> = get_test_r1cs(power);
    let mut rng = ark_std::test_rng();
    let ccs = CCS::<Fr>::from_r1cs(r1cs);
    let (pedersen_params, _) = Pedersen::<Projective>::setup(&mut rng, ccs.n - ccs.l - 1).unwrap();
    let z_1 = get_test_z(power);
    let z_2 = get_test_z(power);

    let (running_instance, w1) = ccs
        .to_lcccs::<_, _, Pedersen<Projective>>(&mut rng, &pedersen_params, &z_1)
        .unwrap();

    let (new_instance, w2) = ccs
        .to_cccs::<_, _, Pedersen<Projective>>(&mut rng, &pedersen_params, &z_2)
        .unwrap();

    let poseidon_config = poseidon_canonical_config::<Fr>();

    let mut transcript_p: PoseidonSponge<Fr> = PoseidonSponge::<Fr>::new(&poseidon_config);

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

    // let mut transcript_v: PoseidonSponge<Fr> = PoseidonSponge::<Fr>::new(&poseidon_config);
    // transcript_v.absorb(&Fr::from_le_bytes_mod_order(b"init init"));
    //
    // // Run the verifier side of the multifolding
    // let folded_lcccs_v = NIMFS::<Projective, PoseidonSponge<Fr>>::verify(
    //     &mut transcript_v,
    //     &ccs,
    //     &[running_instance.clone()],
    //     &[new_instance.clone()],
    //     proof,
    // )
    //     .unwrap();
    // assert_eq!(folded_lcccs, folded_lcccs_v);
    //
    // // Check that the folded LCCCS instance is a valid instance with respect to the folded witness
    // folded_lcccs.check_relation(&ccs, &folded_witness).unwrap();
}

fn main() {
    // let pows: Vec<usize> = (10..24).collect();
    let pows: Vec<usize> = vec![16, 20];
    let iter = 10;
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
    if let Err(e) = write_to_csv(&pows, &prove_times, format!("hypernova_prove_times.csv")) {
        eprintln!("Failed to write to CSV: {}", e);
    } else {
        println!("CSV file has been successfully written.");
    }
}
