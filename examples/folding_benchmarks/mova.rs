use ark_pallas::{Fr, Projective};
use ark_std::log2;
use ark_std::UniformRand;
use folding_schemes::commitment::pedersen::Pedersen;
use folding_schemes::commitment::CommitmentScheme;
use folding_schemes::folding::mova::nifs::NIFS;
use folding_schemes::folding::mova::Witness;
use folding_schemes::transcript::poseidon::poseidon_canonical_config;
use std::mem::size_of_val;
use std::time::{Duration, Instant};

use crate::bench_utils::{get_test_r1cs, get_test_z, write_to_csv};
use folding_schemes::folding::mova::traits::MovaR1CS;

use ark_crypto_primitives::sponge::poseidon::PoseidonSponge;
use ark_crypto_primitives::sponge::CryptographicSponge;
use folding_schemes::arith::r1cs::R1CS;

mod bench_utils;

#[allow(non_snake_case)]
fn mova_benchmark(power: usize, prove_times: &mut Vec<Duration>) {
    let size = 1 << power;

    // define r1cs and parameters
    let r1cs: R1CS<Fr> = get_test_r1cs(power);
    let mut rng = ark_std::test_rng();
    let (pedersen_params, _) = Pedersen::<Projective>::setup(&mut rng, r1cs.A.n_cols).unwrap();
    let poseidon_config = poseidon_canonical_config::<ark_pallas::Fr>();
    let mut transcript_p: PoseidonSponge<Fr> = PoseidonSponge::<Fr>::new(&poseidon_config);

    // INSTANCE 1
    let z_1: Vec<Fr> = get_test_z(power);

    let (w_1, x_1) = r1cs.split_z(&z_1);

    let mut witness_1 = Witness::<Projective>::new(w_1.clone(), r1cs.A.n_rows);
    let vector = vec![1; size];

    // Populate error vector
    witness_1.E = vector.into_iter().map(Fr::from).collect();

    // generate a random evaluation point for MLE
    let size_rE_1 = log2(size);
    let rE_1: Vec<_> = (0..size_rE_1).map(|_| Fr::rand(&mut rng)).collect();

    let committed_instance_1 = witness_1
        .commit::<Pedersen<Projective>>(&pedersen_params, x_1, rE_1)
        .unwrap();

    // INSTANCE 2
    let z_2 = get_test_z(power);
    let (w_2, x_2) = r1cs.split_z(&z_2);
    let mut witness_2 = Witness::<Projective>::new(w_2.clone(), r1cs.A.n_rows);

    let vector = vec![2; size];

    // Populate error vector
    witness_2.E = vector.into_iter().map(Fr::from).collect();

    let size_rE_2 = log2(size);
    let rE_2: Vec<_> = (0..size_rE_2).map(|_| Fr::rand(&mut rng)).collect();

    let committed_instance_2 = witness_2
        .commit::<Pedersen<Projective>>(&pedersen_params, x_2, rE_2)
        .unwrap();

    let start = Instant::now();
    // NIFS.P
    let result = NIFS::<Projective, Pedersen<Projective>, PoseidonSponge<Fr>>::prove(
        &pedersen_params,
        &r1cs,
        &mut transcript_p,
        &committed_instance_1,
        &committed_instance_2,
        &witness_1,
        &witness_2,
    )
    .unwrap();

    let prove_time = start.elapsed();
    prove_times.push(prove_time);
    println!("Mova prove time {:?}", prove_time);
    println!("Mova bytes used {:?}", size_of_val(&result));

    //NIFS.V
    let mut transcript_v: PoseidonSponge<Fr> = PoseidonSponge::<Fr>::new(&poseidon_config);

    let (proof, instance_witness) = result;

    let folded_committed_instance =
        NIFS::<Projective, Pedersen<Projective>, PoseidonSponge<Fr>>::verify(
            &mut transcript_v,
            &committed_instance_1,
            &committed_instance_2,
            &proof,
        )
        .unwrap();
    let check =
        r1cs.check_relaxed_instance_relation(&instance_witness.w, &folded_committed_instance);
    match check {
        Ok(_) => println!("The relation check was successful."),
        Err(e) => println!("The relation check failed: {:?}", e),
    }
}

fn main() {
    // let pows: Vec<usize> = (10..24).collect();
    let pows: Vec<usize> = vec![16, 20];
    let iter = 1;
    let mut prove_times: Vec<Duration> = Vec::with_capacity(pows.len() * iter);
    for i in 0..iter {
        println!("starting {:}", i);

        println!("{:?}", pows);

        for pow in &pows {
            println!("{}", pow);
            mova_benchmark(*pow, &mut prove_times);
        }

        println!("Powers {:?}", pows);
        println!("Prove times {:?}", prove_times);
    }
    if let Err(e) = write_to_csv(&pows, &prove_times, "mova_prove_times.csv".to_string()) {
        eprintln!("Failed to write to CSV: {}", e);
    } else {
        println!("CSV file has been successfully written.");
    }
}
