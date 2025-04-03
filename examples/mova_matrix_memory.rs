use ark_crypto_primitives::sponge::poseidon::PoseidonSponge;
use ark_crypto_primitives::sponge::CryptographicSponge;
use ark_pallas::{Fr, Projective};
use ark_std::{log2, UniformRand};
use folding_schemes::commitment::hyrax::HyraxGenerators;
use folding_schemes::commitment::pedersen::Pedersen;
use folding_schemes::commitment::{CommitmentScheme, NethermindCommitmentScheme};
use folding_schemes::folding::nova::nifs::mova_matrix::{RelaxedCommittedRelation, Witness, NIFS};
use folding_schemes::transcript::poseidon::poseidon_canonical_config;
use folding_schemes::Curve;
use matrex::Matrix;
use rand::{Rng, RngCore};
use std::alloc::GlobalAlloc;
use std::time::{Duration, Instant};

const NUM_OF_PRECONDITION_FOLDS: &[usize] = &[1];

fn random_sparse_matrix<C: Curve>(n: usize, rng: &mut impl RngCore) -> Matrix<C::ScalarField> {
    let elements = (0..n)
        .map(|row| {
            (
                row * n + rand::thread_rng().gen_range(0..n),
                C::ScalarField::rand(rng),
            )
        })
        .collect();
    Matrix::sparse_from_vec(elements, n, n).unwrap()
}

// Helper functions
fn get_instances<C: Curve>(
    num: usize,
    n: usize,
    rng: &mut impl RngCore,
    params: &HyraxGenerators<C>,
) -> Vec<(Witness<C>, RelaxedCommittedRelation<C>)> {
    let start = Instant::now();
    (0..num)
        .map(|_| -> (Witness<C>, RelaxedCommittedRelation<C>) {
            let before_gen_matrices = start.elapsed();
            println!("before_gen_matrices 2 {:?}", before_gen_matrices);
            // A matrix
            let a = random_sparse_matrix::<C>(n, rng);
            // B matrix
            let b = random_sparse_matrix::<C>(n, rng);
            // C = A * B matrix
            let c = (&a * &b).unwrap();
            // Error matrix initialized to 0s
            let e = Matrix::zero(n, n);
            let after_gen_matrices = start.elapsed();
            println!("after_gen_matrices 2 {:?}", after_gen_matrices);

            // Random challenge
            let rE = (0..2 * log2(n))
                .map(|_| C::ScalarField::rand(rng))
                .collect();
            // Witness
            let witness = Witness::new::<false>(a, b, c, e);
            let before_commit_wit = start.elapsed();
            println!("before_commit_wit 2 {:?}", before_commit_wit);
            let instance = witness.commit(params, rE).unwrap();
            let after_commit_wit = start.elapsed();
            println!("after_commit_wit 2 {:?}", after_commit_wit);
            (witness, instance)
        })
        .collect()
}

fn bench_mova_matrix() {
    let mut rng = ark_std::test_rng();
    let mat_dim = 1 << 15; // 4x4 matrices
    println!("mat_dim {}", mat_dim);
    for count in NUM_OF_PRECONDITION_FOLDS {
        println!("Starting with pedersen setup");

        let start = Instant::now();
        let hyrax_params =
            HyraxGenerators::<Projective>::setup(&mut rng, log2(mat_dim * mat_dim) as usize);
        let hyrax_elapsed = start.elapsed();
        println!("hyrax_elapsed 1 {:?}", hyrax_elapsed);

        let poseidon_config = poseidon_canonical_config::<Fr>();
        let pp_hash = Fr::rand(&mut rng);

        let mut total_duration = Duration::ZERO;
        println!("Starting with gen instances");

        let before_instances = start.elapsed();
        println!("before_instances 1 {:?}", before_instances);
        let mut instances: Vec<(Witness<Projective>, RelaxedCommittedRelation<Projective>)> =
            get_instances::<Projective>(
                count + 1, // we want the number of folds plus one for the acc_instance
                mat_dim,
                &mut rng,
                &hyrax_params,
            );
        let after_instances = start.elapsed();
        println!("after_instances 1 {:?}", after_instances);

        let mut transcript_p = PoseidonSponge::<Fr>::new(&poseidon_config);
        let mut acc = instances.pop().unwrap();

        for _ in 0..*count {
            let mut next = instances.pop().unwrap();
            total_duration += {
                let timer = Instant::now();
                println!("Starting with prove");

                let before_proof = start.elapsed();
                println!("before_proof 1 {:?}", before_proof);
                let (wit_acc, inst_acc, _) = NIFS::<Projective, PoseidonSponge<Fr>>::prove(
                    &mut transcript_p,
                    pp_hash,
                    &mut next.0,
                    &next.1,
                    &acc.0,
                    &acc.1,
                )
                .unwrap();
                let after_proof = start.elapsed();
                println!("after_proof 1 {:?}", after_proof);
                let time = timer.elapsed();
                acc = (wit_acc, inst_acc);
                time
            };
        }
        println!("Ending");
    }
}

fn main() {
    bench_mova_matrix();
}
