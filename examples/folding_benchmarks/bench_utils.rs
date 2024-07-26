use ark_ff::PrimeField;
use csv::Writer;
use folding_schemes::arith::r1cs::R1CS;
use folding_schemes::utils::vec::SparseMatrix;
use num_bigint::BigUint;
use rand::Rng;
use std::env;
use std::error::Error;
use std::time::Duration;

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
        coeffs,
    }
}

pub fn get_test_r1cs<F: PrimeField>(power: usize) -> R1CS<F> {
    // Define matrices A, B, and C as specified
    let a = create_large_diagonal_matrix::<F>(power);
    let b = create_large_diagonal_matrix::<F>(power);
    let c = create_large_diagonal_matrix::<F>(power);

    // Return the R1CS structure
    R1CS {
        l: 1,
        A: a,
        B: b,
        C: c,
    }
}

pub fn get_test_z<F: PrimeField>(power: usize) -> Vec<F> {
    let z_vec = create_random_biguints(1 << power, 384);
    to_f_vec(z_vec)
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

pub fn to_f_vec<F: PrimeField>(z: Vec<BigUint>) -> Vec<F> {
    let mut result = Vec::with_capacity(z.len());
    for bigint in z {
        result.push(F::from(bigint));
    }
    result
}

pub fn write_to_csv(
    pows: &[usize],
    prove_times: &[Duration],
    file_path: String,
) -> Result<(), Box<dyn Error>> {
    let path = env::current_dir()?
        .join("examples/folding_benchmarks")
        .join(file_path);
    let mut writer = Writer::from_path(path)?;

    writer.write_record(["pow", "prove_time"])?;

    let mut pows_cycle = pows.iter().cycle();

    for prove_time in prove_times {
        if let Some(pow) = pows_cycle.next() {
            writer.write_record(&[pow.to_string(), prove_time.as_micros().to_string()])?;
        }
    }

    writer.flush()?;
    Ok(())
}
