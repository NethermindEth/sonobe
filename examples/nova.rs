use ark_ff::PrimeField;
use ark_pallas::{Fr, Projective};
use ark_std::UniformRand;
use folding_schemes::ccs::r1cs::R1CS;
use folding_schemes::commitment::pedersen::Pedersen;
use folding_schemes::commitment::CommitmentScheme;
use folding_schemes::folding::nova::nifs::NIFS;
use folding_schemes::folding::nova::traits::NovaR1CS;
use folding_schemes::folding::nova::Witness;
use folding_schemes::utils::vec::{dense_matrix_to_sparse, SparseMatrix};
use rand::Rng;
use std::mem::size_of_val;
use std::time::Instant;

fn main() {
    println!("starting");
    let random_num: usize = rand::thread_rng().gen_range(1..=2642245);
    let r1cs = get_test_r1cs();
    let z = get_test_z(random_num);
    let (w, x) = r1cs.split_z(&z);

    let mut rng = ark_std::test_rng();
    let (pedersen_params, _) = Pedersen::<Projective>::setup(&mut rng, r1cs.A.n_cols).unwrap();

    let running_instance_w = Witness::<Projective>::new(w.clone(), r1cs.A.n_rows);
    let running_committed_instance = running_instance_w
        .commit::<Pedersen<Projective>>(&pedersen_params, x)
        .unwrap();

    let incoming_instance_z = get_test_z(4);
    let (w, x) = r1cs.split_z(&incoming_instance_z);
    let incoming_instance_w = Witness::<Projective>::new(w.clone(), r1cs.A.n_rows);
    let incoming_committed_instance = incoming_instance_w
        .commit::<Pedersen<Projective>>(&pedersen_params, x)
        .unwrap();

    let r = Fr::rand(&mut rng); // folding challenge would come from the RO

    let start = Instant::now();
    // NIFS.P
    let (T, cmT) = NIFS::<Projective, Pedersen<Projective>>::compute_cmT(
        &pedersen_params,
        &r1cs,
        &running_instance_w,
        &running_committed_instance,
        &incoming_instance_w,
        &incoming_committed_instance,
    )
    .unwrap();
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

    println!("Z number: {:?}", random_num);
    println!("Nova prove time {:?}", start.elapsed());
    println!("Nova bytes used {:?}", size_of_val(&result));

    let (folded_w, _) = result;

    let folded_committed_instance = NIFS::<Projective, Pedersen<Projective>>::verify(
        r,
        &running_committed_instance,
        &incoming_committed_instance,
        &cmT,
    );
    r1cs.check_relaxed_instance_relation(&folded_w, &folded_committed_instance)
        .unwrap();
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
