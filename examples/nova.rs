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
use std::time::Instant;

mod bench_utils;

fn main() {
    println!("starting");

    let big_num = BigUint::one() << 100;
    let r1cs = get_test_r1cs();

    let z = get_test_z_albert(big_num);

    let (w, x) = r1cs.split_z(&z);

    let mut rng = ark_std::test_rng();
    let (pedersen_params, _) = Pedersen::<Projective>::setup(&mut rng, r1cs.A.n_cols).unwrap();

    let running_instance_w = Witness::<Projective>::new(w.clone(), r1cs.A.n_rows);
    let running_committed_instance = running_instance_w
        .commit::<Pedersen<Projective>>(&pedersen_params, x)
        .unwrap();

    let four = BigUint::one() + BigUint::one() + BigUint::one() + BigUint::one();
    let incoming_instance_z = get_test_z_albert(four);
    let (w, x) = r1cs.split_z(&incoming_instance_z);
    let incoming_instance_w = Witness::<Projective>::new(w.clone(), r1cs.A.n_rows);
    let incoming_committed_instance = incoming_instance_w
        .commit::<Pedersen<Projective>>(&pedersen_params, x)
        .unwrap();

    let poseidon_config = poseidon_canonical_config::<ark_pallas::Fr>();
    let mut transcript_p = PoseidonTranscript::<Projective>::new(&poseidon_config);
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

    println!("Nova prove time {:?}", start.elapsed());
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
