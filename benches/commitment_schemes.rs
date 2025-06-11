use ark_ec::VariableBaseMSM;
use ark_pallas::{Fr, Projective};
use ark_std::{UniformRand, Zero};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use folding_schemes::commitment::hyrax::{Hyrax, HyraxGenerators};
use folding_schemes::commitment::pedersen::Pedersen;
use folding_schemes::commitment::{CommitmentScheme, NethermindCommitmentScheme};
use folding_schemes::Curve;
use matrex::Matrix;
use rand::{Rng, RngCore};
use std::time::Duration;

// Matrix sizes to test (dimensions will be 2^size x 2^size)
const MATRIX_SIZES: &[usize] = &[4, 6, 8, 10]; // 16, 64, 128, 1024
const SPARSITY_FACTOR: f64 = 0.05; // 5% non-zero elements for sparse matrices

/// Generates a random dense matrix of specified dimensions
fn random_dense_matrix<C: Curve>(n: usize, rng: &mut impl RngCore) -> Matrix<C::ScalarField> {
    let elements: Vec<C::ScalarField> = (0..n * n).map(|_| C::ScalarField::rand(rng)).collect();
    Matrix::dense_from_vec(elements, n, n).unwrap()
}

/// Generates a random sparse matrix with specified dimensions and sparsity
fn random_sparse_matrix<C: Curve>(n: usize, rng: &mut impl RngCore) -> Matrix<C::ScalarField> {
    let non_zero_count = (n * n * SPARSITY_FACTOR as usize).max(1);

    let mut elements = Vec::with_capacity(non_zero_count);
    for _ in 0..non_zero_count {
        let row = rand::thread_rng().gen_range(0..n);
        let col = rand::thread_rng().gen_range(0..n);
        let index = row * n + col;
        let value = C::ScalarField::rand(rng);
        elements.push((index, value));
    }

    Matrix::sparse_from_vec(elements, n, n).unwrap()
}

/// Benchmark Pedersen vs Hyrax for dense matrices
fn bench_dense_commits(c: &mut Criterion) {
    let mut group = c.benchmark_group("Dense Matrix Commitments");
    group.measurement_time(Duration::from_secs(10));

    let mut rng = ark_std::test_rng();

    for &size in MATRIX_SIZES {
        let n = 1 << size; // 2^size
        let matrix = random_dense_matrix::<Projective>(n, &mut rng);
        let data = matrix.as_dense_slice().unwrap();

        // Setup parameters
        let pedersen_params = Pedersen::<Projective>::setup_prover(&mut rng, n * n).unwrap();
        let hyrax_params = HyraxGenerators::<Projective>::setup(&mut rng, n * n);

        group.bench_with_input(BenchmarkId::new("Pedersen", n), &n, |b, _| {
            b.iter(|| Pedersen::<Projective>::commit(&pedersen_params, data, &Fr::zero()).unwrap());
        });

        group.bench_with_input(BenchmarkId::new("Hyrax", n), &n, |b, _| {
            b.iter(|| Hyrax::<Projective>::commit(data, &hyrax_params).unwrap());
        });
    }

    group.finish();
}

/// Benchmark Pedersen vs Hyrax for sparse matrices
fn bench_sparse_commits(c: &mut Criterion) {
    let mut group = c.benchmark_group("Sparse Matrix Commitments");
    group.measurement_time(Duration::from_secs(10));

    let mut rng = ark_std::test_rng();

    for &size in MATRIX_SIZES {
        let n = 1 << size; // 2^size
        let matrix = random_sparse_matrix::<Projective>(n, &mut rng);
        let sparse_data = matrix.as_sparse_slice().unwrap();

        // Setup parameters
        let pedersen_params = Pedersen::<Projective>::setup_prover(&mut rng, n * n).unwrap();
        let hyrax_params = HyraxGenerators::<Projective>::setup(&mut rng, n * n);

        group.bench_with_input(BenchmarkId::new("Pedersen Sparse", n), &n, |b, _| {
            b.iter(|| {
                Pedersen::<Projective>::commit_sparse(&pedersen_params, sparse_data, &Fr::zero())
                    .unwrap()
            });
        });

        group.bench_with_input(BenchmarkId::new("Hyrax Sparse", n), &n, |b, _| {
            b.iter(|| {
                Hyrax::<Projective>::commit_sparse_matrix(sparse_data, &hyrax_params).unwrap()
            });
        });
    }

    group.finish();
}
/// Benches the same MSM in 2 scenarios: a single MSM, or multiple smaller MSM until reaching the same size
fn bench_msm(c: &mut Criterion) {
    let mut group = c.benchmark_group("MultiScalar Exponentiations");
    group.measurement_time(Duration::from_secs(10));
    let mut rng = ark_std::test_rng();

    for &size in MATRIX_SIZES {
        let n = 1 << size; // 2^size
        let matrix = random_dense_matrix::<Projective>(n, &mut rng);
        let data = matrix.as_dense_slice().unwrap();
        let partial_data: Vec<_> = (0..n).map(|_| Fr::rand(&mut rng)).collect();

        // Use fully qualified type syntax for Affine
        type ProjectiveAffine = <Projective as ark_ec::CurveGroup>::Affine;

        // Setup parameters
        let generators: Vec<ProjectiveAffine> =
            (0..n).map(|_| ProjectiveAffine::rand(&mut rng)).collect();
        let partial_gen = &generators[..n];

        let pedersen_params = Pedersen::<Projective>::setup_prover(&mut rng, n * n).unwrap();

        group.bench_with_input(BenchmarkId::new("Single Large MSM", n), &n, |b, _| {
            b.iter(|| Projective::msm_unchecked(&pedersen_params.generators, data));
        });

        group.bench_with_input(BenchmarkId::new("Many Small MSM", n), &n, |b, _| {
            b.iter(|| {
                // Do n MSM os size n instead of a large MSM osf size n*n.
                for _ in 0..n {
                    let _ = Projective::msm_unchecked(&partial_gen, &partial_data);
                }
            });
        });
    }
    group.finish();
}

criterion_group!(
    benches,
    bench_msm,
    bench_dense_commits,
    bench_sparse_commits
);
criterion_main!(benches);
