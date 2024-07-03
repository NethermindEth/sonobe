use ark_ff::{PrimeField, BigInteger};
use num_bigint::BigUint;
use num_traits::One;
use folding_schemes::ccs::r1cs::R1CS;
use folding_schemes::Error;
use folding_schemes::utils::vec::{dense_matrix_to_sparse, SparseMatrix};

// would be best to move this to other file
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

pub fn get_test_z<F: PrimeField>(input: BigUint) -> Vec<F> {
    let one = BigUint::one();
    let five = &one + &one + &one + &one + &one;

    to_F_vec(vec![
        one,                                  // 1
        input.clone(),                             // io
        &input * &input * &input + &input + &five, // x^3 + x + 5
        &input * &input,                           // x^2
        &input * &input * &input,                  // x^3
        &input * &input * &input + &input,         // x^3 + x
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
pub fn to_F_vec<F: PrimeField>(z: Vec<BigUint>) -> Vec<F> {
    let mut result = Vec::with_capacity(z.len()); // Pre-allocate space for efficiency
    for bigint in z {
        // Convert each BigUint to F::BigInt
        // match F::try_from(bigint) {
        match num_bigint_to_ark_bigint::<F>(&bigint) {
            Ok(f_bigint) => {
                // Attempt to convert F::BigInt to the prime field element
                if let Some(field_element) = F::from_bigint(f_bigint) {
                    result.push(field_element);
                } else {
                    // Handle the case where the conversion is not possible
                    eprintln!("Conversion to field element failed for bigint: {:?}", f_bigint);
                    continue; // Optionally skip or handle differently
                }
                // result.push(f_bigint);
            },
            Err(e) => {
                // Handle errors from bigint conversion
                eprintln!("Error converting bigint: {:?}", e);
                continue; // Optionally skip or handle differently
            }
        }
    }
    result
}

pub fn num_bigint_to_ark_bigint<F: PrimeField>(value: &BigUint) -> Result<F::BigInt, Error> {
    F::BigInt::try_from(value.clone()).map_err(|_| {
        Error::BigIntConversionError("Failed to convert to PrimeField::BigInt".to_string())
    })
}

// pub const BIG_NUM: BigUint = BigUint::one() << 80;
// pub const FOUR: BigUint = BigUint::one() + BigUint::one() + BigUint::one() + BigUint::one();

