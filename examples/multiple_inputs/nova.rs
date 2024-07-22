use crate::bench_utils::{get_test_r1cs, get_test_z_albert};
use ark_ff::{ BigInteger, Field, PrimeField};
use ark_pallas::{Fr, Projective};
use ark_std::{log2, UniformRand};
use folding_schemes::commitment::pedersen::Pedersen;
use folding_schemes::commitment::CommitmentScheme;
use folding_schemes::folding::nova::nifs::NIFS;
use folding_schemes::folding::nova::traits::NovaR1CS;
use folding_schemes::folding::nova::Witness;
use folding_schemes::transcript::poseidon::{poseidon_canonical_config, PoseidonTranscript};
use folding_schemes::transcript::Transcript;
use folding_schemes::utils::sum_check::{ SumCheck};
use num_bigint::BigUint;
use num_traits::{One, Zero};
use rand::Rng;
use std::mem::size_of_val;
use std::time::{Duration, Instant};

use std::error::Error;
use csv::Writer;

mod bench_utils;

fn nova_benchmark(power: usize, prove_times: &mut Vec<Duration>) {
    let big_num: BigUint = BigUint::one() << 100;
    let r1cs = get_test_r1cs(power);

    let z = get_test_z_albert(big_num.clone(), power);

    let (w, x) = r1cs.split_z(&z);

    let mut rng = ark_std::test_rng();
    let (pedersen_params, _) = Pedersen::<Projective>::setup(&mut rng, r1cs.A.n_cols).unwrap();

    let running_instance_w = Witness::<Projective>::new(w.clone(), r1cs.A.n_rows);
    let running_committed_instance = running_instance_w
        .commit::<Pedersen<Projective>>(&pedersen_params, x)
        .unwrap();

    let four = BigUint::one() + BigUint::one() + BigUint::one() + BigUint::one();
    let incoming_instance_z = get_test_z_albert(four, power);
    let (w, x) = r1cs.split_z(&incoming_instance_z);
    let incoming_instance_w = Witness::<Projective>::new(w.clone(), r1cs.A.n_rows);
    let incoming_committed_instance = incoming_instance_w
        .commit::<Pedersen<Projective>>(&pedersen_params, x)
        .unwrap();

    let poseidon_config = poseidon_canonical_config::<ark_pallas::Fr>();
    let mut transcript_p = PoseidonTranscript::<Projective>::new(&poseidon_config);
    // let vector = vec![1; size];
    // //
    // witness_1.E = vector.into_iter().map(|x| Fr::from(x)).collect();
    //
    // let vector = vec![2; size];
    // //
    // witness_2.E = vector.into_iter().map(|x| Fr::from(x)).collect();
    // NIFS.P
    let start = Instant::now();

    let (T, cmT) = NIFS::<Projective, Pedersen<Projective>>::compute_cmT(
        &pedersen_params,
        &r1cs,
        &running_instance_w,
        &running_committed_instance,
        &incoming_instance_w,
        &incoming_committed_instance,
    )
    .unwrap();

    match transcript_p.absorb_point(&cmT) {
        Ok(_) => {
            //
        }
        Err(e) => {
            println!("Absorbed failed: {:?}", e);
        }
    }

    let r = transcript_p.get_challenge();
    let result = NIFS::<Projective, Pedersen<Projective>>::fold_instances(
        r,
        &running_instance_w,
        &running_committed_instance,
        &incoming_instance_w,
        &incoming_committed_instance,
        &T,
        cmT,
    )
    .unwrap();

    let prove_time = start.elapsed();
    prove_times.push(prove_time);
    println!("Nova prove time {:?}", prove_time);
    println!("Nova bytes used {:?}", size_of_val(&result));

    let (folded_w, _) = result;

    let folded_committed_instance = NIFS::<Projective, Pedersen<Projective>>::verify(
        r,
        &running_committed_instance,
        &incoming_committed_instance,
        &cmT,
    );
    let check = r1cs.check_relaxed_instance_relation(&folded_w, &folded_committed_instance);
    match check {
        Ok(_) => println!("The relation check was successful."),
        Err(e) => println!("The relation check failed: {:?}", e),
    }
}

fn write_to_csv(prove_times: &mut Vec<Duration>) -> Result<(), Box<dyn Error>> {
    let mut wtr = Writer::from_path("prove_times.csv")?;
    wtr.write_record(&["x", "y", "z"])?;
    wtr.flush()?;
    Ok(())
}


fn main() {
    println!("starting");

    let pows: Vec<usize> = (10..24).collect();
    println!("{:?}", pows);

    let mut prove_times: Vec<Duration> = Vec::new();

    for pow in pows {
        nova_benchmark(pow, &mut prove_times);
    }

    let pows_print: Vec<usize> = (10..24).collect();
    println!("Powers {:?}", pows_print);

    println!("Prove times {:?}", prove_times);

    println!(
        "| {0: <10} | {1: <10} |",
        "2^pow", "prove time"
    );
    for i in 0..prove_times.len() {
        println!("| {0: <10} | {1:?} |", pows_print[i], prove_times[i]);
    }

    // write_to_csv(&mut prove_times);
    
}
