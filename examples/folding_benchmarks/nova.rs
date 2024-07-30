use crate::bench_utils::{get_test_r1cs, get_test_z, write_to_csv};
use ark_pallas::{Fr, Projective};
use folding_schemes::commitment::pedersen::Pedersen;
use folding_schemes::commitment::CommitmentScheme;
use folding_schemes::folding::nova::nifs::NIFS;
use folding_schemes::folding::nova::traits::NovaR1CS;
use folding_schemes::folding::nova::{Witness};
use folding_schemes::transcript::poseidon::poseidon_canonical_config;
use folding_schemes::transcript::Transcript;
use std::time::{Duration, Instant};

use ark_crypto_primitives::sponge::poseidon::PoseidonSponge;
use ark_crypto_primitives::sponge::CryptographicSponge;
use ark_std::UniformRand;

mod bench_utils;

fn nova_benchmark(power: usize, prove_times: &mut Vec<Duration>) {
    let size = 1 << power;

    let r1cs = get_test_r1cs(power);

    let z_1 = get_test_z(power, 384);


    let (w, x) = r1cs.split_z(&z_1);

    let mut rng = ark_std::test_rng();
    let (pedersen_params, _) = Pedersen::<Projective>::setup(&mut rng, r1cs.A.n_cols).unwrap();

    let mut witness_1 = Witness::<Projective>::new(w.clone(), r1cs.A.n_rows);
    let vector = vec![1; size];
    //
    witness_1.E = vector.into_iter().map(Fr::from).collect();
    let mut running_committed_instance = witness_1
        .commit::<Pedersen<Projective>>(&pedersen_params, x)
        .unwrap();

    running_committed_instance.u = Fr::rand(&mut rng);

    let z_2 = get_test_z(power, 20);

    let (w, x) = r1cs.split_z(&z_2);
    let mut witness_2 = Witness::<Projective>::new(w.clone(), r1cs.A.n_rows);
    let vector = vec![2; size];
    //
    witness_2.E = vector.into_iter().map(Fr::from).collect();
    let mut incoming_committed_instance = witness_2
        .commit::<Pedersen<Projective>>(&pedersen_params, x)
        .unwrap();

    incoming_committed_instance.u = Fr::rand(&mut rng);

    let poseidon_config = poseidon_canonical_config::<Fr>();
    let mut transcript_p: PoseidonSponge<Fr> = PoseidonSponge::<Fr>::new(&poseidon_config);



    // NIFS.P
    let start = Instant::now();

    let (t, cm_t) = NIFS::<Projective, Pedersen<Projective>>::compute_cmT(
        &pedersen_params,
        &r1cs,
        &witness_1,
        &running_committed_instance,
        &witness_2,
        &incoming_committed_instance,
    )
        .unwrap();

    let elapsed = start.elapsed();
    println!("Time before Randomness generation {:?}", elapsed);
    transcript_p.absorb_nonnative(&cm_t);

    let r = transcript_p.get_challenge();
    let elapsed = start.elapsed();
    println!("Time aftre Randomness generation {:?}", elapsed);

    let elapsed = start.elapsed();
    println!("Time before starting folding {:?}", elapsed);

    let result = NIFS::<Projective, Pedersen<Projective>>::fold_instances(
        r,
        &witness_1,
        &running_committed_instance,
        &witness_2,
        &incoming_committed_instance,
        &t,
        cm_t,
    )
        .unwrap();
    let elapsed = start.elapsed();

    println!("Time after folding {:?}", elapsed);

    let prove_time = start.elapsed();
    prove_times.push(prove_time);
    println!("Nova prove time {:?}", prove_time);

    let (folded_w, _) = result;

    let folded_committed_instance = NIFS::<Projective, Pedersen<Projective>>::verify(
        r,
        &running_committed_instance,
        &incoming_committed_instance,
        &cm_t,
    );
    let check = r1cs.check_relaxed_instance_relation(&folded_w, &folded_committed_instance);
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
            nova_benchmark(*pow, &mut prove_times);
        }

        println!("Powers {:?}", pows);
        println!("Prove times {:?}", prove_times);
    }
    if let Err(e) = write_to_csv(&pows, &prove_times, "nova_prove_times.csv".to_string()) {
        eprintln!("Failed to write to CSV: {}", e);
    } else {
        println!("CSV file has been successfully written.");
    }
}
