use ark_ff::{BigInteger, PrimeField};
use folding_schemes::ccs::r1cs::R1CS;
use folding_schemes::utils::vec::{dense_matrix_to_sparse, SparseMatrix};
use folding_schemes::Error;
use num_bigint::BigUint;
use num_traits::One;
use rand::Rng;

// would be best to move this to other file
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
const POWER: usize = 16;
fn create_large_diagonal_matrix<F: PrimeField>() -> SparseMatrix<F> {
    let size = 1 << POWER;
    let mut coeffs: Vec<Vec<(F, usize)>> = Vec::with_capacity(size);

    // Populate the diagonal elements
    for i in 0..size {
        // Each row has one non-zero entry at (i, i) with a value of 2
        coeffs.push(vec![(F::from(1u64), i)]);
    }

    // Instantiate SparseMatrix directly
    SparseMatrix {
        n_rows: size,
        n_cols: size,
        coeffs
    }
}

pub fn get_test_r1cs<F: PrimeField>() -> R1CS<F> {
    // Define matrices A, B, and C as specified
    let A = create_large_diagonal_matrix::<F>();
    let B = create_large_diagonal_matrix::<F>();
    let C = create_large_diagonal_matrix::<F>();

    // let A = create_large_diagonal_matrix_2::<F>();
    // let B = create_large_diagonal_matrix_2::<F>();
    // let C = create_large_diagonal_matrix_2::<F>();

    // println!("A: {:?} {:?}", A.coeffs, A.n_rows);
    // println!("B: {:?} {:?}", B.coeffs, B.n_rows);
    // println!("C: {:?} {:?}", C.coeffs, C.n_rows);

    // Return the R1CS structure
    R1CS { l: 1, A, B, C }
}

// pub fn get_test_z<F: PrimeField>(input: BigUint) -> Vec<F> {
//     let one = BigUint::one();
//     let five = &one + &one + &one + &one + &one;
//     to_F_vec_2(vec![
//         one,                                       // 1
//         input.clone(),                             // io
//         &input * &input * &input + &input + &five, // x^3 + x + 5
//         &input * &input,                           // x^2
//         &input * &input * &input,                  // x^3
//         &input * &input * &input + &input,
//         &input * &input * &input + &input, // x^3 + x
//         &input * &input * &input + &input, // x^3 + x
//                                            // x^3 + x
//     ])
// }

pub fn get_test_z<F: PrimeField>(input: BigUint) -> Vec<F> {
    let one = BigUint::one();
    let five = &one + &one + &one + &one + &one;
    let size = 1 << POWER; // Calculate size only once
    let mut z_vec = Vec::with_capacity(size); // Preallocate memory for efficiency

    // Add initial elements
    z_vec.push(one); // 1
    z_vec.push(input.clone()); // input

    // Pre-compute input cubed since it does not change in the loop
    let input_cubed = &input * &input * &input;

    // Fill the rest of the vector
    z_vec.extend(std::iter::repeat(input_cubed).take(size - 2));
    to_F_vec_2(z_vec)
}

pub fn get_test_z_albert<F: PrimeField>(input: BigUint) -> Vec<F> {
    // let one = BigUint::one();
    // let five = &one + &one + &one + &one + &one;
    // let size = 1 << POWER; // Calculate size only once
    // let mut z_vec = Vec::with_capacity(size); // Preallocate memory for // efficiency
//
    // // Add initial elements
    // // z_vec.push(one); // 1
    // z_vec.push(input.clone()); // input
//
    // // Pre-compute input cubed since it does not change in the loop
    // let input_cubed = &input * &input * &input;
//
    // // Fill the rest of the vector
    // // z_vec.extend(std::iter::repeat(input_cubed).take(size - 2));
    // z_vec.extend((0..size - 1).map(|i| &input_cubed * (i as u32) + i*i*i*i + i*i + i*i*i));


    let z_vec = create_random_biguints(1<<POWER, 384);
    // for (i, num) in z_vec.iter().enumerate() {
    //     // println!("Random BigUint {}: {}", i + 1, num);
    // }
    to_F_vec_2(z_vec)
}

pub fn create_random_biguints(count: usize, max_bits: u32) -> Vec<BigUint> {
    let mut rng = rand::thread_rng();

    (0..count)
        .map(|_| {
            let bit_size = rng.gen_range(1..=max_bits);
            let mut bytes = vec![0u8; (bit_size as usize + 7) / 8];
            rng.fill(&mut bytes[..]);
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
                    eprintln!(
                        "Conversion to field element failed for bigint: {:?}",
                        f_bigint
                    );
                    continue; // Optionally skip or handle differently
                }
                // result.push(f_bigint);
            }
            Err(e) => {
                // Handle errors from bigint conversion
                eprintln!("Error converting bigint: {:?}", e);
                continue; // Optionally skip or handle differently
            }
        }
    }
    result
}
pub fn to_F_vec_2<F: PrimeField>(z: Vec<BigUint>) -> Vec<F> {
    let mut result = Vec::with_capacity(z.len()); // Pre-allocate space for efficiency
    for bigint in z {
        // Convert each BigUint to F::BigInt
        match F::try_from(bigint) {
            // match num_bigint_to_ark_bigint::<F>(&bigint) {
            Ok(f_bigint) => {
                // // Attempt to convert F::BigInt to the prime field element
                // if let Some(field_element) = F::from_bigint(f_bigint) {
                //     result.push(field_element);
                // } else {
                //     // Handle the case where the conversion is not possible
                //     eprintln!("Conversion to field element failed for bigint: {:?}", f_bigint);
                //     continue; // Optionally skip or handle differently
                // }
                result.push(f_bigint);
            }
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

fn create_large_diagonal_matrix_2<F: PrimeField>() -> SparseMatrix<F> {
    let size = 1 << POWER; // 2^16
    let mut matrix = vec![vec![0; size]; size];

    for i in 0..size {
        matrix[i][i] = 2; // Set diagonal elements to 2
    }

    // println!("{:?}", matrix);
    to_F_matrix(matrix)
}
//
// pub fn get_test_r1cs<F: PrimeField>() -> R1CS<F> {
//     // Define matrices A, B, and C as specified
//     let A = create_large_diagonal_matrix::<F>();
//     let B = create_large_diagonal_matrix::<F>();
//     let C = create_large_diagonal_matrix::<F>();
//     println!("A: {:?} {:?}", A.n_cols, A.n_rows);
//     println!("B: {:?} {:?}", B.n_cols, B.n_rows);
//     println!("C: {:?} {:?}", C.n_cols, C.n_rows);
//
//     // Return the R1CS structure
//     R1CS::<F> { l: 1, A, B, C }
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