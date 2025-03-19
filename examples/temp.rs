use ark_crypto_primitives::sponge::poseidon::PoseidonSponge;
use ark_crypto_primitives::sponge::CryptographicSponge;
use ark_ec::CurveGroup;
use ark_pallas::{Fr, Projective};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{log2, UniformRand};
use criterion::{criterion_group, criterion_main, Criterion};
use folding_schemes::commitment::pedersen::Pedersen;
use folding_schemes::commitment::{CommitmentScheme, SparseCommitmentScheme};
use folding_schemes::folding::nova::nifs::mova_matrix::{RelaxedCommittedRelation, Witness, NIFS};
use folding_schemes::transcript::poseidon::poseidon_canonical_config;
use folding_schemes::{Curve, Error};
use matrex::{Matrix, SparseMatrix};
use rand::{Rng, RngCore};
use sha3::digest::{ExtendableOutput, Update};
use sha3::Shake256;
use std::io::Read;
use std::time::{Duration, Instant};
use num_integer::Roots;
use rand_chacha::ChaCha20Rng;
use rand_chacha::rand_core::SeedableRng;
use folding_schemes::utils::mle::MultilinearExtension;


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

#[derive(Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct PedersenGenerators<G: CurveGroup> {
    pub generators: Vec<G>,
}

impl<G: CurveGroup> PedersenGenerators<G> {
    #[tracing::instrument(skip_all, name = "PedersenGenerators::new")]
    pub fn new(len: usize, label: &[u8]) -> Self {
        let mut shake = Shake256::default();
        shake.update(label);
        let mut buf = vec![];
        G::generator().serialize_compressed(&mut buf).unwrap();
        shake.update(&buf);

        let mut reader = shake.finalize_xof();
        let mut seed = [0u8; 32];
        reader.read_exact(&mut seed).unwrap();
        let mut rng = ChaCha20Rng::from_seed(seed);

        let mut generators: Vec<G> = Vec::new();
        for _ in 0..len {
            generators.push(G::rand(&mut rng));
        }

        Self { generators }
    }
}

// Helper functions
fn get_instances<C: Curve, CS: CommitmentScheme<C> + SparseCommitmentScheme<C>>(
    num: usize,
    n: usize,
    rng: &mut impl RngCore,
    params: &CS::ProverParams,
) -> Vec<(Witness<C>, RelaxedCommittedRelation<C>)> {
    (0..num)
        .map(|_| -> (Witness<C>, RelaxedCommittedRelation<C>) {
            // A matrix
            let a = random_sparse_matrix::<C>(n, rng);
            println!("a {:?}", a.len());

            // B matrix
            let b = random_sparse_matrix::<C>(n, rng);
            // C = A * B matrix
            let c = (a.clone() * &b).unwrap();
            // Error matrix initialized to 0s
            let e = Matrix::zero(n, n);
            // Random challenge
            let rE = (0..2 * log2(n))
                .map(|_| C::ScalarField::rand(rng))
                .collect();
            // Witness
            let witness = Witness::new::<false>(a, b, c, e);
            let instance = witness.commit::<CS, false>(params, rE).unwrap();
            (witness, instance)
        })
        .collect()
}

fn bench_mova_matrix() {
    let mut rng = ark_std::test_rng();
    let mat_dim = 1<<14; // 4x4 matrices
    println!("{}", mat_dim);
    let a1 = random_sparse_matrix::<Projective>(mat_dim , &mut rng);
    // B matrix
    let b1 = random_sparse_matrix::<Projective>(mat_dim , &mut rng);

    let a2 = random_sparse_matrix::<Projective>(mat_dim , &mut rng);
    // B matrix
    let b2 = random_sparse_matrix::<Projective>(mat_dim , &mut rng);
    let c1 = random_sparse_matrix::<Projective>(mat_dim , &mut rng);
    // B matrix
    let c2 = random_sparse_matrix::<Projective>(mat_dim , &mut rng);
    let u = Fr::rand(&mut  rng);


    let A1B2 = (&a1* &b2).unwrap();

    let B1A2 = (&a2 * &b1).unwrap();
    let A1B2B1A2 = (A1B2 + B1A2).unwrap(); // Sparse (but less sparse)

    let u2c1 = c1.clone() * u;
    let T = ((A1B2B1A2 - &c2).unwrap() - u2c1).unwrap();
    let rE: Vec<_> = (0..2 * log2((mat_dim * mat_dim).sqrt()))
        .map(|_| Fr::rand(&mut rng))
        .collect();
    // println!("T {:?}", T);
    let mle = MultilinearExtension::from_evaluations(&T, log2(a1.len()) as usize);
    let mleT_evaluated = mle.evaluate(&rE);




    // for count in NUM_OF_PRECONDITION_FOLDS {
    //     // Set up transcript and commitment scheme
    //     println!("Starting with pedersen setup");
    //     let pedersen_params =
    //         Pedersen::<Projective>::setup2(&mut rng, mat_dim * mat_dim).unwrap();
    //     // println!("Pedersen::setup() called!");
    //     // println!("{:?}", std::backtrace::Backtrace::capture());
    //     // println!("Starting jolt setup");
    //
    //     // let temp2 = PedersenGenerators::<Projective>::new(mat_dim * mat_dim, b"nick");
    //     let poseidon_config = poseidon_canonical_config::<Fr>();
    //     let pp_hash = Fr::rand(&mut rng);
    //
    //     let mut total_duration = Duration::ZERO;
    //     println!("Starting with gen instances");
    //
    //     let mut instances: Vec<(Witness<Projective>, RelaxedCommittedRelation<Projective>)> =
    //         get_instances::<Projective, Pedersen<Projective>>(
    //             count + 1, // we want the number of folds plus one for the acc_instance
    //             mat_dim,
    //             &mut rng,
    //             &pedersen_params,
    //         );
    //     let mut transcript_p = PoseidonSponge::<Fr>::new(&poseidon_config);
    //     let mut acc = instances.pop().unwrap();
    //
    //     for _ in 0..*count {
    //         let next = instances.pop().unwrap();
    //         total_duration += {
    //             let timer = Instant::now();
    //             println!("Starting with prove");
    //
    //             let (wit_acc, inst_acc, _) =
    //                 NIFS::<Projective, Pedersen<Projective>, PoseidonSponge<Fr>>::prove(
    //                     &mut transcript_p,
    //                     pp_hash,
    //                     &next.0,
    //                     &next.1,
    //                     &acc.0,
    //                     &acc.1,
    //                 )
    //                 .unwrap();
    //             let time = timer.elapsed();
    //             acc = (wit_acc, inst_acc);
    //             time
    //         };
    //     }
    // }
}

fn main() {
    bench_mova_matrix();

}
