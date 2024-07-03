use std::marker::PhantomData;
use ark_ff::{Field, PrimeField};
use ark_pallas::{Fr, Projective};
use ark_std::{log2, UniformRand};
use folding_schemes::ccs::r1cs::R1CS;
use folding_schemes::commitment::pedersen::Pedersen;
use folding_schemes::commitment::CommitmentScheme;
use folding_schemes::folding::nova::nifs::NIFS;
use folding_schemes::folding::nova::traits::NovaR1CS;
use folding_schemes::folding::nova::Witness;
use folding_schemes::transcript::poseidon::{poseidon_canonical_config, PoseidonTranscript};
use folding_schemes::transcript::Transcript;
use folding_schemes::utils::vec::{dense_matrix_to_sparse, SparseMatrix};
use rand::Rng;
use std::mem::size_of_val;
use std::sync::Arc;
use std::time::Instant;
use num_bigint::BigUint;
use crate::bench_utils::{get_test_r1cs, get_test_z};
use num_traits::{One, Zero};
use folding_schemes::Error;
use folding_schemes::utils::mle::dense_vec_to_dense_mle;
use folding_schemes::utils::sum_check::{IOPSumCheck, SumCheck};
use folding_schemes::utils::virtual_polynomial::{build_eq_x_r_vec, VirtualPolynomial, VPAuxInfo};


mod bench_utils;


fn main() {
    println!("starting");

    let big_num = BigUint::one() << 80;
    let r1cs = get_test_r1cs();

    let z = get_test_z(big_num);

    let (w, x) = r1cs.split_z(&z);


    let mut rng = ark_std::test_rng();
    let (pedersen_params, _) = Pedersen::<Projective>::setup(&mut rng, r1cs.A.n_cols).unwrap();

    let running_instance_w = Witness::<Projective>::new(w.clone(), r1cs.A.n_rows);
    let running_committed_instance = running_instance_w
        .commit::<Pedersen<Projective>>(&pedersen_params, x)
        .unwrap();

    let four = BigUint::one() + BigUint::one() + BigUint::one() + BigUint::one();
    let incoming_instance_z = get_test_z(four);
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

    // let g = compute_g(running_committed_instance, incoming_committed_instance, running_instance_w, incoming_instance_w, &transcript_p.get_challenge())?;

    // let vars = log2(running_instance_w.E.len()) as usize;
    //
    // let mut g: VirtualPolynomial<Fr> = VirtualPolynomial::new(2);
    // let temp = [Fr::zero()];
    // let mle = dense_vec_to_dense_mle(2, &temp);
    //
    // let vp_aux_info = VPAuxInfo::<Fr> {
    //     max_degree: 2,
    //     num_variables: 2,
    //     phantom: PhantomData::<Fr>,
    // };
    //
    // g.add_mle_list([Arc::new(mle.clone()), Arc::new(mle.clone())], Fr::one()).expect("TODO: panic message");
    //
    // let sumcheck_proof = IOPSumCheck::<Projective,
    //     PoseidonTranscript<Projective>>::prove(&g, &mut transcript_p)
    //     .map_err(|err| Error::SumCheckProveError(err.to_string())).unwrap();
    //
    // let sumcheck_subclaim =
    //     IOPSumCheck::<Projective, PoseidonTranscript<Projective>>::verify(running_committed_instance.x[0], &sumcheck_proof, &vp_aux_info, &mut transcript_p)
    //         .map_err(|err| Error::SumCheckVerifyError(err.to_string()));


}

// pub fn get_test_r1cs<F: PrimeField>() -> R1CS<F> {
//     // R1CS for: x^3 + x + 5 = y (example from article
//     // https://www.vitalik.ca/general/2016/12/10/qap.html )
//     let A = to_F_matrix::<F>(vec![
//         vec![0, 1, 0, 0, 0, 0],
//         vec![0, 0, 0, 1, 0, 0],
//         vec![0, 1, 0, 0, 1, 0],
//         vec![5, 0, 0, 0, 0, 1],
//     ]);
//     let B = to_F_matrix::<F>(vec![
//         vec![0, 1, 0, 0, 0, 0],
//         vec![0, 1, 0, 0, 0, 0],
//         vec![1, 0, 0, 0, 0, 0],
//         vec![1, 0, 0, 0, 0, 0],
//     ]);
//     let C = to_F_matrix::<F>(vec![
//         vec![0, 0, 0, 1, 0, 0],
//         vec![0, 0, 0, 0, 1, 0],
//         vec![0, 0, 0, 0, 0, 1],
//         vec![0, 0, 1, 0, 0, 0],
//     ]);
//
//     R1CS::<F> { l: 1, A, B, C }
// }
//
// pub fn get_test_z<F: PrimeField>(input: usize) -> Vec<F> {
//     // z = (1, io, w)
//     to_F_vec(vec![
//         1,
//         input,                             // io
//         input * input * input + input + 5, // x^3 + x + 5
//         input * input,                     // x^2
//         input * input * input,             // x^2 * x
//         input * input * input + input,     // x^3 + x
//     ])
// }
//
// pub fn to_F_matrix<F: PrimeField>(M: Vec<Vec<usize>>) -> SparseMatrix<F> {
//     dense_matrix_to_sparse(to_F_dense_matrix(M))
// }
// pub fn to_F_dense_matrix<F: PrimeField>(M: Vec<Vec<usize>>) -> Vec<Vec<F>> {
//     M.iter()
//         .map(|m| m.iter().map(|r| F::from(*r as u64)).collect())
//         .collect()
// }
// pub fn to_F_vec<F: PrimeField>(z: Vec<usize>) -> Vec<F> {
//     z.iter().map(|c| F::from(*c as u64)).collect()
// }
