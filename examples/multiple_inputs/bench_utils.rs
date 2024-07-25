use std::env;
use std::error::Error;
use std::time::Duration;
use ark_ff::{ PrimeField};
use csv::Writer;
use folding_schemes::arith::r1cs;
use folding_schemes::utils::vec::{dense_matrix_to_sparse, SparseMatrix};
use num_bigint::BigUint;
use rand::Rng;
use folding_schemes::arith::r1cs::R1CS;

fn create_large_diagonal_matrix<F: PrimeField>(power: usize) -> SparseMatrix<F> {
    let size = 1 << power;
    let mut coeffs: Vec<Vec<(F, usize)>> = Vec::with_capacity(size);

    // Populate the diagonal elements
    for i in 0..size {
        // Each row has one non-zero entry at (i, i) with a value of 1
        coeffs.push(vec![(F::from(1u64), i)]);
    }
    // Instantiate SparseMatrix directly
    SparseMatrix {
        n_rows: size,
        n_cols: size,
        coeffs
    }
}

pub fn get_test_r1cs<F: PrimeField>(power: usize) -> R1CS<F> {
    // Define matrices A, B, and C as specified
    let A = create_large_diagonal_matrix::<F>(power);
    let B = create_large_diagonal_matrix::<F>(power);
    let C = create_large_diagonal_matrix::<F>(power);

    // Return the R1CS structure
    R1CS { l: 1, A, B, C }
}

pub fn get_test_z<F: PrimeField>(power: usize) -> Vec<F> {

    let z_vec = create_random_biguints(1<< power, 384);
    to_F_vec(z_vec)
}

pub fn create_random_biguints(size: usize, max_bits: u32) -> Vec<BigUint> {
    let mut rng = rand::thread_rng();
    let mut bytes = vec![0u8; (max_bits as usize + 7) / 8];

    (0..size)
        .map(|_| {
            // Generate a random bit size between 1 and max_bits
            let bit_size = rng.gen_range(1..=max_bits);
            // Calculate the number of bytes needed for the given bit size
            let byte_size = (bit_size as usize + 7) / 8;

            // Fill the bytes vector with random data up to the byte_size
            rng.fill(&mut bytes[..byte_size]);

            // Convert the byte slice to a BigUint and return it
            BigUint::from_bytes_le(&bytes) % (BigUint::from(1u32) << bit_size)
        })
        .collect()
}

pub fn to_F_matrix<F: PrimeField>(M: Vec<Vec<usize>>) -> SparseMatrix<F> {
    dense_matrix_to_sparse(to_F_dense_matrix(M))
}
pub fn to_F_dense_matrix<F: PrimeField>(M: Vec<Vec<usize>>) -> Vec<Vec<F>> {
    M.iter()
        .map(|m| m.iter().map(|r| F::from(*r as u64)).collect())
        .collect()
}
pub fn to_F_vec<F: PrimeField>(z: Vec<BigUint>) -> Vec<F> {
    let mut result = Vec::with_capacity(z.len());
    for bigint in z {
        // Convert each BigUint to F::BigInt
        match F::try_from(bigint) {
            Ok(f_bigint) => result.push(f_bigint),
            Err(e) => eprintln!("Error converting bigint: {:?}", e)
        }
    }
    result
}

pub fn write_to_csv(pows: &[usize], prove_times: &[Duration], file_path: String) -> Result<(), Box<dyn Error>> {
    let path = env::current_dir()?.join("examples/multiple_inputs").join(file_path);
    let mut writer = Writer::from_path(path)?;

    writer.write_record(&["pow", "prove_time"])?;

    for (pow, prove_time) in pows.iter().zip(prove_times) {
        writer.write_record(&[
            pow.to_string(),
            prove_time.as_micros().to_string(),
        ])?;
    }

    writer.flush()?;
    Ok(())
}