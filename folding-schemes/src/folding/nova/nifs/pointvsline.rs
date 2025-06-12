use crate::folding::nova::nifs::mova::{CommittedInstance, Witness};
use crate::folding::nova::nifs::mova_matrix::{RelaxedCommittedRelation, Witness as MatrixWitness};
use crate::transcript::Transcript;
use crate::utils::mle::{dense_vec_to_dense_mle, MultilinearExtension, SparseOrDensePolynomial};
use crate::{Curve, Error};
use ark_crypto_primitives::sponge::Absorb;
use ark_ff::{One, PrimeField};
use ark_poly::univariate::{DensePolynomial, SparsePolynomial};
use ark_poly::{DenseMultilinearExtension, DenseUVPolynomial, Polynomial};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{log2, Zero};
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};
use std::fmt::Debug;

/// Implements the Points vs Line as described in
/// [Mova](https://eprint.iacr.org/2024/1220.pdf) and Section 4.5.2 from Thaler’s book

pub struct PointVsLineEvaluationClaimR1CS<C: Curve> {
    pub mleE1_prime: C::ScalarField,
    pub mleE2_prime: C::ScalarField,
    pub rE_prime: Vec<C::ScalarField>,
}
/// Proof from step 1 protocol 6
#[derive(Debug, Clone, Eq, PartialEq, CanonicalSerialize, CanonicalDeserialize)]
pub struct PointVsLineProofR1CS<C: Curve> {
    pub h1: DensePolynomial<C::ScalarField>,
    pub h2: DensePolynomial<C::ScalarField>,
}

pub struct PointVsLineEvaluationClaimMatrix<C: Curve> {
    pub mleE2_prime: C::ScalarField,
    pub rE_prime: Vec<C::ScalarField>,
}
/// Proof from step 1 protocol 6
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct PointVsLineProofMatrix<C: Curve> {
    pub h2: SparseOrDensePolynomial<C::ScalarField>,
}

pub trait PointVsLine<C: Curve, T: Transcript<C::ScalarField>> {
    type PointVsLineProof: Debug + Clone;

    type PointVsLineEvaluationClaim;

    type CommittedInstance: Debug + Clone + Absorb; // + CommittedInstanceOps<C>;
    type Witness: Debug + Clone;

    fn prove(
        transcript: &mut T,
        ci1: Option<&Self::CommittedInstance>,
        ci2: &Self::CommittedInstance,
        w1: &Self::Witness,
        w2: &Self::Witness,
    ) -> Result<(Self::PointVsLineProof, Self::PointVsLineEvaluationClaim), Error>;

    fn verify(
        transcript: &mut T,
        ci1: Option<&Self::CommittedInstance>,
        ci2: &Self::CommittedInstance,
        proof: &Self::PointVsLineProof,
        mleE1_prime: Option<&C::ScalarField>,
        mleE2_prime: &<C>::ScalarField,
        rE_prime_p: &[<C>::ScalarField], // the rE_prime of the prover
    ) -> Result<
        Vec<<C>::ScalarField>, // rE=rE1'=rE2'.
        Error,
    >;
}
#[derive(Clone, Debug, Default)]
pub struct PointVsLineR1CS<C: Curve, T: Transcript<C::ScalarField>> {
    _phantom_C: std::marker::PhantomData<C>,
    _phantom_T: std::marker::PhantomData<T>,
}

/// Protocol 6 from Mova
impl<C: Curve, T: Transcript<C::ScalarField>> PointVsLine<C, T> for PointVsLineR1CS<C, T> {
    type PointVsLineProof = PointVsLineProofR1CS<C>;
    type PointVsLineEvaluationClaim = PointVsLineEvaluationClaimR1CS<C>;
    type CommittedInstance = CommittedInstance<C>;
    type Witness = Witness<C>;

    fn prove(
        transcript: &mut T,
        ci1: Option<&Self::CommittedInstance>,
        ci2: &Self::CommittedInstance,
        w1: &Self::Witness,
        w2: &Self::Witness,
    ) -> Result<(Self::PointVsLineProof, Self::PointVsLineEvaluationClaim), Error> {
        let ci1 = ci1.ok_or_else(|| Error::Other("Missing ci1 in R1CS prove".to_string()))?;

        let n_vars: usize = log2(w1.E.len()) as usize;

        let mleE1 = dense_vec_to_dense_mle(n_vars, &w1.E);
        let mleE2 = dense_vec_to_dense_mle(n_vars, &w2.E);

        // We have l(0) = r1, l(1) = r2 so we know that l(x) = r1 + x(r2-r1) thats why we need r2-r1
        let r2_sub_r1: Vec<<C>::ScalarField> = ci1
            .rE
            .iter()
            .zip(&ci2.rE)
            .map(|(&r1, r2)| *r2 - r1)
            .collect();

        let h1 = compute_h(&mleE1, &ci1.rE, &r2_sub_r1)?;
        let h2 = compute_h(&mleE2, &ci1.rE, &r2_sub_r1)?;

        transcript.absorb(&h1.coeffs());
        transcript.absorb(&h2.coeffs());

        let beta = transcript.get_challenge();

        let mleE1_prime = h1.evaluate(&beta);
        let mleE2_prime = h2.evaluate(&beta);

        let rE_prime = compute_l(&ci1.rE, &r2_sub_r1, beta)?;

        Ok((
            Self::PointVsLineProof { h1, h2 },
            Self::PointVsLineEvaluationClaim {
                mleE1_prime,
                mleE2_prime,
                rE_prime,
            },
        ))
    }

    fn verify(
        transcript: &mut T,
        ci1: Option<&Self::CommittedInstance>,
        ci2: &Self::CommittedInstance,
        proof: &Self::PointVsLineProof,
        mleE1_prime: Option<&C::ScalarField>,
        mleE2_prime: &<C>::ScalarField,
        rE_prime_p: &[<C>::ScalarField],
    ) -> Result<Vec<<C>::ScalarField>, Error> {
        let ci1 = ci1.ok_or_else(|| Error::Other("Missing ci1 in R1CS verify".to_string()))?;
        if proof.h1.evaluate(&C::ScalarField::zero()) != ci1.mleE {
            return Err(Error::NotEqual);
        }

        if proof.h2.evaluate(&C::ScalarField::one()) != ci2.mleE {
            return Err(Error::NotEqual);
        }

        transcript.absorb(&proof.h1.coeffs());
        transcript.absorb(&proof.h2.coeffs());

        let beta = transcript.get_challenge();

        if let Some(mleE1_prime_val) = mleE1_prime {
            if *mleE1_prime_val != proof.h1.evaluate(&beta) {
                return Err(Error::NotEqual);
            }
        } else {
            return Err(Error::Other(
                "Missing mleE1_prime in R1CS verify".to_string(),
            ));
        }

        if *mleE2_prime != proof.h2.evaluate(&beta) {
            return Err(Error::NotEqual);
        }

        let r2_sub_r1: Vec<<C>::ScalarField> = ci1
            .rE
            .iter()
            .zip(&ci2.rE)
            .map(|(&r1, r2)| *r2 - r1)
            .collect();
        let rE_prime = compute_l(&ci1.rE, &r2_sub_r1, beta)?;
        if rE_prime != rE_prime_p {
            return Err(Error::NotEqual);
        }

        Ok(rE_prime)
    }
}

#[derive(Clone, Debug, Default)]
pub struct PointVsLineMatrix<C: Curve, T: Transcript<C::ScalarField>> {
    _phantom_C: std::marker::PhantomData<C>,
    _phantom_T: std::marker::PhantomData<T>,
}

impl<C: Curve, T: Transcript<C::ScalarField>> PointVsLine<C, T> for PointVsLineMatrix<C, T> {
    type PointVsLineProof = PointVsLineProofMatrix<C>;
    type PointVsLineEvaluationClaim = PointVsLineEvaluationClaimMatrix<C>;
    type CommittedInstance = RelaxedCommittedRelation<C>;
    type Witness = MatrixWitness<C>;

    fn prove(
        transcript: &mut T,
        _ci1: Option<&Self::CommittedInstance>,
        ci2: &Self::CommittedInstance,
        _w1: &Self::Witness,
        w2: &Self::Witness,
    ) -> Result<(Self::PointVsLineProof, Self::PointVsLineEvaluationClaim), Error> {
        // Derive randomness
        let r1_scalar = C::ScalarField::from_le_bytes_mod_order(b"r1");
        transcript.absorb(&r1_scalar);
        let r1: Vec<C::ScalarField> = transcript.get_challenges(ci2.rE.len());

        let n_vars: usize = log2(w2.E.len()) as usize;
        let mleE2 = MultilinearExtension::from_evaluations(&w2.E, n_vars);

        // We have l(0) = r1, l(1) = r2 so we know that l(x) = r1 + x(r2-r1) that's why we need r2-r1
        let r2_sub_r1: Vec<<C>::ScalarField> =
            r1.iter().zip(&ci2.rE).map(|(&r1, r2)| *r2 - r1).collect();

        let h2 = compute_h2(&mleE2, &r1, &r2_sub_r1)?;

        transcript.absorb(&h2.coeffs());

        let beta = transcript.get_challenge();

        let mleE2_prime = h2.evaluate(&beta);

        let rE_prime = compute_l(&r1, &r2_sub_r1, beta)?;

        Ok((
            Self::PointVsLineProof { h2 },
            Self::PointVsLineEvaluationClaim {
                mleE2_prime,
                rE_prime,
            },
        ))
    }

    fn verify(
        transcript: &mut T,
        _ci1: Option<&Self::CommittedInstance>,
        ci2: &Self::CommittedInstance,
        proof: &Self::PointVsLineProof,
        _mleE1_prime: Option<&C::ScalarField>,
        mleE2_prime: &<C>::ScalarField,
        rE_prime_p: &[<C>::ScalarField],
    ) -> Result<Vec<<C>::ScalarField>, Error> {
        if proof.h2.evaluate(&C::ScalarField::one()) != ci2.mleE {
            return Err(Error::NotEqual);
        }

        let r1_scalar = C::ScalarField::from_le_bytes_mod_order(b"r1");
        transcript.absorb(&r1_scalar);

        let r1 = transcript.get_challenges(ci2.rE.len());

        transcript.absorb(&proof.h2.coeffs());

        let beta = transcript.get_challenge();

        if *mleE2_prime != proof.h2.evaluate(&beta) {
            return Err(Error::NotEqual);
        }

        let r2_sub_r1: Vec<<C>::ScalarField> =
            r1.iter().zip(&ci2.rE).map(|(&r1, r2)| *r2 - r1).collect();
        let rE_prime = compute_l(&r1, &r2_sub_r1, beta)?;
        if rE_prime != rE_prime_p {
            return Err(Error::NotEqual);
        }

        Ok(rE_prime)
    }
}

fn compute_h<F: PrimeField>(
    mle: &DenseMultilinearExtension<F>,
    r1: &[F],
    r2_sub_r1: &[F],
) -> Result<DensePolynomial<F>, Error> {
    let n_vars: usize = mle.num_vars;
    if r1.len() != r2_sub_r1.len() || r1.len() != n_vars {
        return Err(Error::NotEqual);
    }

    // Start with coefficient vectors. For now they are constant polynomials with a single coefficient
    let mut coeffs: Vec<Vec<F>> = mle.evaluations.iter().map(|&x| vec![x]).collect();

    for (i, (&r1_i, &r2_sub_r1_i)) in r1.iter().zip(r2_sub_r1.iter()).enumerate().take(n_vars) {
        // Create a linear polynomial r(X) = r1_i + (r2_sub_r1_i) * X (basically l)
        let half_len = 1 << (n_vars - i - 1);
        let new_coeffs: Vec<Vec<F>> = (0..half_len)
            .into_par_iter()
            .map(|b| {
                let left_idx = b << 1;
                let right_idx = left_idx + 1;

                let left_coeffs = &coeffs[left_idx];
                let right_coeffs = &coeffs[right_idx];

                // Initialize result coefficients
                let mut result_coeffs = vec![F::zero(); right_coeffs.len() + 1];

                // Add left polynomial contribution
                for (j, &left_val) in left_coeffs.iter().enumerate() {
                    result_coeffs[j] = left_val;
                }

                // Add (right - left) * (r1_i + r2_sub_r1_i * X) contribution directly
                for (j, (&right_val, &left_val)) in
                    right_coeffs.iter().zip(left_coeffs.iter()).enumerate()
                {
                    let diff = right_val - left_val;
                    result_coeffs[j] += diff * r1_i;
                    result_coeffs[j + 1] += diff * r2_sub_r1_i;
                }

                result_coeffs
            })
            .collect();

        coeffs = new_coeffs;
    }

    Ok(DensePolynomial::from_coefficients_vec(
        coeffs.swap_remove(0),
    ))
}

/// Implementation for computing h by not following Algorithm 1 "MLE-after-line composition" off the Mova paper
/// This is due to the need to support sparse representation.
/// Currently, this is only used for the mova_matrix.rs implementation configured to use Matrex
fn compute_h2<F: PrimeField>(
    mle: &MultilinearExtension<F>,
    r1: &[F],
    r2_sub_r1: &[F],
) -> Result<SparseOrDensePolynomial<F>, Error> {
    let n_vars = mle.num_vars();

    if r1.len() != r2_sub_r1.len() || r1.len() != n_vars {
        return Err(Error::NotEqual);
    }

    match mle {
        MultilinearExtension::DenseMLE(mle_dense) => {
            // Start with evaluations as degree-0 constant polynomials,
            // We'll represent polynomials as coefficient vectors instead of DensePolynomials as it's more efficient.
            let mut coeffs: Vec<Vec<F>> = mle_dense
                .evaluations
                .iter()
                .map(|&eval| vec![eval])
                .collect();

            // For each variable, fold pairs of polynomials
            for (i, (&r1_i, &r2_sub_r1_i)) in r1.iter().zip(r2_sub_r1.iter()).enumerate() {
                let half_len = 1 << (n_vars - i - 1);

                let new_coeffs: Vec<Vec<F>> = (0..half_len)
                    .into_par_iter()
                    .map(|b| {
                        let left_idx = b << 1;
                        let right_idx = left_idx + 1;

                        let left_coeffs: &Vec<F> = &coeffs[left_idx];
                        let right_coeffs: &Vec<F> = &coeffs[right_idx];

                        let max_degree = right_coeffs.len() + 1;
                        let mut result_coeffs = vec![F::zero(); max_degree];

                        // Add left polynomial first
                        for (j, &left_val) in left_coeffs.iter().enumerate() {
                            result_coeffs[j] = left_val;
                        }

                        // Add right polynomial contribution directly: (right - left) * (r1_i + r2_sub_r1_i * x)
                        for (j, (&right_val, &left_val)) in
                            right_coeffs.iter().zip(left_coeffs.iter()).enumerate()
                        {
                            let diff = right_val - left_val;
                            result_coeffs[j] += diff * r1_i;
                            result_coeffs[j + 1] += diff * r2_sub_r1_i;
                        }

                        result_coeffs
                    })
                    .collect();

                coeffs = new_coeffs;
            }

            // Convert final coefficient vector to polynomial
            Ok(SparseOrDensePolynomial::from_dense(
                DensePolynomial::from_coefficients_vec(coeffs.into_iter().next().unwrap()),
            ))
        }

        MultilinearExtension::SparseMLE(mle_sparse) => {
            // If there are no evaluations, return the zero polynomial
            if mle_sparse.evaluations.is_empty() {
                return Ok(SparseOrDensePolynomial::from_sparse(
                    SparsePolynomial::zero(),
                ));
            }
            let max_degree = n_vars + 1;
            // Pre-compute linear factors to avoid repeated computation
            let linear_factors: Vec<(F, F, F, F)> = (0..n_vars)
                .map(|i| {
                    (
                        r1[i],            // factor_1_const
                        r2_sub_r1[i],     // factor_1_linear
                        F::one() - r1[i], // factor_0_const
                        -r2_sub_r1[i],    // factor_0_linear
                    )
                })
                .collect();

            let result_coeffs = mle_sparse
                .evaluations
                .par_iter()
                .map(|(&index, &value)| {
                    let mut contrib_coeffs = vec![F::zero(); max_degree];
                    contrib_coeffs[0] = value;
                    let mut current_degree = 0;

                    // Multiply by the linear factor for each variable
                    for i in 0..n_vars {
                        let bit_i = (index >> i) & 1 == 1;
                        let (const_term, linear_term) = if bit_i {
                            // If bit_i == 1, use r1_i + r2_sub_r1_i * x
                            (linear_factors[i].0, linear_factors[i].1)
                        } else {
                            // If bit_i == 0, use 1 - r1_i - r2_sub_r1_i * x
                            (linear_factors[i].2, linear_factors[i].3)
                        };

                        // Multiply in-place by linear polynomial
                        contrib_coeffs[current_degree + 1] =
                            contrib_coeffs[current_degree] * linear_term;
                        for j in (1..=current_degree).rev() {
                            contrib_coeffs[j] = contrib_coeffs[j] * const_term
                                + contrib_coeffs[j - 1] * linear_term;
                        }
                        contrib_coeffs[0] *= const_term;

                        current_degree += 1;
                    }

                    // Return just the required coefficients
                    contrib_coeffs.truncate(current_degree + 1);
                    contrib_coeffs
                })
                .reduce(
                    || vec![F::zero(); max_degree],
                    |mut acc, contrib| {
                        // Parallel reduction: combine two coefficient vectors
                        for (i, &coeff) in contrib.iter().enumerate() {
                            if i < acc.len() {
                                acc[i] += coeff;
                            }
                        }
                        acc
                    },
                );

            // Remove trailing zeros
            let mut result_coeffs = result_coeffs;
            while result_coeffs.len() > 1 && result_coeffs.last() == Some(&F::zero()) {
                result_coeffs.pop();
            }

            Ok(SparseOrDensePolynomial::from_dense(
                DensePolynomial::from_coefficients_vec(result_coeffs),
            ))
        }
    }
}

fn compute_l<F: PrimeField>(r1: &[F], r2_sub_r1: &[F], x: F) -> Result<Vec<F>, Error> {
    if r1.len() != r2_sub_r1.len() {
        return Err(Error::NotEqual);
    }

    // we have l(x) = r1 + x(r2-r1) so return the result
    Ok(r1
        .iter()
        .zip(r2_sub_r1)
        .map(|(&r1, &r2_sub_r1)| r1 + x * r2_sub_r1)
        .collect())
}

#[cfg(test)]
mod tests {
    use super::{
        compute_h, compute_h2, compute_l, PointVsLine, PointVsLineMatrix, PointVsLineR1CS,
    };
    use crate::commitment::pedersen::Pedersen;
    use crate::commitment::CommitmentScheme;
    use crate::transcript::poseidon::poseidon_canonical_config;
    use crate::Error;
    use ark_crypto_primitives::sponge::poseidon::PoseidonSponge;
    use ark_pallas::{Fr, Projective};
    use ark_poly::{
        DenseMultilinearExtension, DenseUVPolynomial, Polynomial, SparseMultilinearExtension,
    };
    use ark_std::{log2, UniformRand};

    use crate::folding::nova::nifs::mova::Witness;
    use crate::folding::nova::nifs::mova_matrix::{
        RelaxedCommittedRelation, Witness as MatrixWitness,
    };

    use crate::commitment::hyrax::{Hyrax, HyraxGenerators};
    use crate::utils::mle::MultilinearExtension;
    use ark_crypto_primitives::sponge::CryptographicSponge;
    use ark_ff::{One, Zero};
    use matrex::Matrix;

    #[test]
    fn test_compute_h() -> Result<(), Error> {
        let mle = DenseMultilinearExtension::from_evaluations_slice(1, &[Fr::from(1), Fr::from(2)]);
        let r0 = [Fr::from(5)];
        let r1 = [Fr::from(6)];
        let r1_sub_r0: Vec<Fr> = r1.iter().zip(&r0).map(|(&x, y)| x - y).collect();

        let result = compute_h(&mle, &r0, &r1_sub_r0)?;
        assert_eq!(
            result,
            DenseUVPolynomial::from_coefficients_slice(&[Fr::from(6), Fr::from(1)])
        );

        let mle = DenseMultilinearExtension::from_evaluations_slice(1, &[Fr::from(1), Fr::from(2)]);
        let r0 = [Fr::from(4)];
        let r1 = [Fr::from(7)];
        let r1_sub_r0: Vec<Fr> = r1.iter().zip(&r0).map(|(&x, y)| x - y).collect();

        let result = compute_h(&mle, &r0, &r1_sub_r0)?;
        assert_eq!(
            result,
            DenseUVPolynomial::from_coefficients_slice(&[Fr::from(5), Fr::from(3)])
        );

        let mle = DenseMultilinearExtension::from_evaluations_slice(
            2,
            &[Fr::from(1), Fr::from(2), Fr::from(3), Fr::from(4)],
        );
        let r0 = [Fr::from(5), Fr::from(4)];
        let r1 = [Fr::from(2), Fr::from(7)];
        let r1_sub_r0: Vec<Fr> = r1.iter().zip(&r0).map(|(&x, y)| x - y).collect();

        let result = compute_h(&mle, &r0, &r1_sub_r0)?;
        assert_eq!(
            result,
            DenseUVPolynomial::from_coefficients_slice(&[Fr::from(14), Fr::from(3)])
        );
        let mle = DenseMultilinearExtension::from_evaluations_slice(
            3,
            &[
                Fr::from(1),
                Fr::from(2),
                Fr::from(3),
                Fr::from(4),
                Fr::from(5),
                Fr::from(6),
                Fr::from(7),
                Fr::from(8),
            ],
        );
        let r0 = [Fr::from(1), Fr::from(2), Fr::from(3)];
        let r1 = [Fr::from(5), Fr::from(6), Fr::from(7)];
        let r1_sub_r0: Vec<Fr> = r1.iter().zip(&r0).map(|(&x, y)| x - y).collect();

        let result = compute_h(&mle, &r0, &r1_sub_r0)?;
        assert_eq!(
            result,
            DenseUVPolynomial::from_coefficients_slice(&[Fr::from(18), Fr::from(28)])
        );
        Ok(())
    }

    #[test]
    fn test_compute_h_errors() {
        let mle = DenseMultilinearExtension::from_evaluations_slice(1, &[Fr::from(1), Fr::from(2)]);
        let r0 = [Fr::from(5)];
        let r1_sub_r0 = [];
        let result = compute_h(&mle, &r0, &r1_sub_r0);
        assert!(result.is_err());

        let mle = DenseMultilinearExtension::from_evaluations_slice(
            2,
            &[Fr::from(1), Fr::from(2), Fr::from(1), Fr::from(2)],
        );
        let r0 = [Fr::from(4)];
        let r1 = [Fr::from(7)];
        let r1_sub_r0: Vec<Fr> = r1.iter().zip(&r0).map(|(&x, y)| x - y).collect();

        let result = compute_h(&mle, &r0, &r1_sub_r0);
        assert!(result.is_err())
    }

    #[test]
    fn test_compute_l() -> Result<(), Error> {
        // Test with simple non-zero values
        let r1 = vec![Fr::from(1), Fr::from(2), Fr::from(3)];
        let r2_sub_r1 = vec![Fr::from(4), Fr::from(5), Fr::from(6)];
        let x = Fr::from(2);

        let expected = vec![
            Fr::from(1) + Fr::from(2) * Fr::from(4),
            Fr::from(2) + Fr::from(2) * Fr::from(5),
            Fr::from(3) + Fr::from(2) * Fr::from(6),
        ];

        let result = compute_l(&r1, &r2_sub_r1, x)?;
        assert_eq!(result, expected);
        Ok(())
    }

    #[test]
    fn test_evaluations_R1CS() -> Result<(), Error> {
        // Basic test with no zero error term to ensure that the folding is correct.
        // This test mainly focuses on if the evaluation of h0 and h1 are correct.
        let mut rng = ark_std::test_rng();

        let (pedersen_params, _) = Pedersen::<Projective>::setup(&mut rng, 4)?;
        let poseidon_config = poseidon_canonical_config::<Fr>();
        let mut transcript_p = PoseidonSponge::<Fr>::new(&poseidon_config);
        let mut transcript_v = PoseidonSponge::<Fr>::new(&poseidon_config);

        let W_i = Witness {
            E: vec![Fr::from(25), Fr::from(50), Fr::from(0), Fr::from(0)],
            W: vec![Fr::from(35), Fr::from(9), Fr::from(27), Fr::from(30)],
            rW: Fr::zero(),
        };
        let rE = (0..log2(W_i.E.len())).map(|_| Fr::rand(&mut rng)).collect();
        // x is not important
        let x = vec![Fr::from(35), Fr::from(9), Fr::from(27), Fr::from(30)];
        let U_i =
            Witness::commit::<Pedersen<Projective>, false>(&W_i, &pedersen_params, x.clone(), rE)?;

        let w_i = Witness {
            E: vec![Fr::from(75), Fr::from(100), Fr::from(0), Fr::from(0)],
            W: vec![Fr::from(35), Fr::from(9), Fr::from(27), Fr::from(30)],
            rW: Fr::zero(),
        };
        let rE = (0..log2(W_i.E.len())).map(|_| Fr::rand(&mut rng)).collect();
        let u_i = Witness::commit::<Pedersen<Projective>, false>(&w_i, &pedersen_params, x, rE)?;

        let (proof, claim) =
            PointVsLineR1CS::prove(&mut transcript_p, Some(&U_i), &u_i, &W_i, &w_i)?;

        let result = PointVsLineR1CS::verify(
            &mut transcript_v,
            Some(&U_i),
            &u_i,
            &proof,
            Some(&claim.mleE1_prime),
            &claim.mleE2_prime,
            &claim.rE_prime,
        );

        assert!(result.is_ok(), "Verification failed");
        // Check if the re_prime is the same
        let re_verified = result.unwrap();
        assert!(re_verified == claim.rE_prime);
        let mut transcript_v = PoseidonSponge::<Fr>::new(&poseidon_config);

        // Pass the wrong committed instance which should result in a wrong evaluation in h returning an error
        let result = PointVsLineR1CS::verify(
            &mut transcript_v,
            Some(&U_i),
            &U_i,
            &proof,
            Some(&claim.mleE1_prime),
            &claim.mleE2_prime,
            &claim.rE_prime,
        );

        assert!(result.is_err(), "Verification was okay when it should fail");

        Ok(())
    }

    #[test]
    fn h2_test_mismatched_input_lengths() {
        let mle = MultilinearExtension::DenseMLE(
            DenseMultilinearExtension::<Fr>::from_evaluations_vec(2, vec![Fr::zero(); 4]),
        );
        let r1 = vec![Fr::one(), Fr::one()];
        let r2_sub_r1 = vec![Fr::one(), Fr::one(), Fr::one()];

        let result = compute_h2(&mle, &r1, &r2_sub_r1);
        assert!(matches!(result, Err(Error::NotEqual)));
    }

    #[test]
    fn test_evaluations_Matrix_dense() -> Result<(), Error> {
        // Basic test with no zero error term to ensure that the folding is correct.
        // This test mainly focuses on if the evaluation of h1 are correct.
        let mut rng = ark_std::test_rng();

        let hyrax_params = HyraxGenerators::<Projective>::setup(&mut rng, 4);
        let poseidon_config = poseidon_canonical_config::<Fr>();
        let mut transcript_p = PoseidonSponge::<Fr>::new(&poseidon_config);
        let mut transcript_v = PoseidonSponge::<Fr>::new(&poseidon_config);

        let three = Fr::one() + Fr::one() + Fr::one();
        let four = three + Fr::one();
        let five = four + Fr::one();
        let six = five + Fr::one();
        let W_i = MatrixWitness {
            A: Matrix::dense_from_vec(vec![three, four, five, six], 2, 2).unwrap(),
            B: Matrix::dense_from_vec(vec![three, four, five, six], 2, 2).unwrap(),
            C: Matrix::dense_from_vec(vec![three, four, five, six], 2, 2).unwrap(),
            E: Matrix::dense_from_vec(vec![three, four, five, six], 2, 2).unwrap(),
        };
        let rE: Vec<Fr> = (0..log2(W_i.E.len())).map(|_| Fr::rand(&mut rng)).collect();
        // below is the commit code modified to work with dense matrices.
        let U_i = {
            let mle = MultilinearExtension::from_evaluations(&W_i.E, log2(W_i.E.len()) as usize);
            let mleE = mle.evaluate(&rE);
            // Right now we are ignoring the hiding property and directly commit to the matrices
            let com_a = Hyrax::commit(W_i.A.as_dense_slice().unwrap(), &hyrax_params)?;
            let com_b = Hyrax::commit(W_i.B.as_dense_slice().unwrap(), &hyrax_params)?;
            let com_c = Hyrax::commit(W_i.C.as_dense_slice().unwrap(), &hyrax_params)?;

            RelaxedCommittedRelation {
                cmA: com_a,
                cmB: com_b,
                cmC: com_c,
                u: Fr::one(),
                mleE,
                rE,
            }
        };

        let w_i = MatrixWitness {
            A: Matrix::dense_from_vec(vec![three, four, five, six], 2, 2).unwrap(),
            B: Matrix::dense_from_vec(vec![three, four, five, six], 2, 2).unwrap(),
            C: Matrix::dense_from_vec(vec![three, four, five, six], 2, 2).unwrap(),
            E: Matrix::dense_from_vec(vec![three, four, five, six], 2, 2).unwrap(),
        };
        let rE: Vec<Fr> = (0..log2(W_i.E.len())).map(|_| Fr::rand(&mut rng)).collect();
        // below is the commit code modified to work with dense matrices.
        let u_i = {
            let mle = MultilinearExtension::from_evaluations(&W_i.E, log2(W_i.E.len()) as usize);
            let mleE = mle.evaluate(&rE);
            // Right now we are ignoring the hiding property and directly commit to the matrices
            let com_a = Hyrax::commit(W_i.A.as_dense_slice().unwrap(), &hyrax_params)?;
            let com_b = Hyrax::commit(W_i.B.as_dense_slice().unwrap(), &hyrax_params)?;
            let com_c = Hyrax::commit(W_i.C.as_dense_slice().unwrap(), &hyrax_params)?;

            RelaxedCommittedRelation {
                cmA: com_a,
                cmB: com_b,
                cmC: com_c,
                u: Fr::one(),
                mleE,
                rE,
            }
        };

        let (proof, claim) = PointVsLineMatrix::prove(&mut transcript_p, None, &u_i, &W_i, &w_i)?;

        let result = PointVsLineMatrix::verify(
            &mut transcript_v,
            None,
            &u_i,
            &proof,
            None,
            &claim.mleE2_prime,
            &claim.rE_prime,
        );

        assert!(result.is_ok(), "Verification failed");
        // Check if the re_prime is the same
        let re_verified = result.unwrap();
        assert!(re_verified == claim.rE_prime);
        let mut transcript_v = PoseidonSponge::<Fr>::new(&poseidon_config);

        // Pass the wrong committed instance which should result in a wrong evaluation in h returning an error
        let result = PointVsLineMatrix::verify(
            &mut transcript_v,
            None,
            &U_i,
            &proof,
            None,
            &claim.mleE2_prime,
            &claim.rE_prime,
        );

        assert!(result.is_err(), "Verification was okay when it should fail");

        Ok(())
    }

    #[test]
    fn test_evaluations_Matrix_sparse() -> Result<(), Error> {
        // Basic test with no zero error term to ensure that the folding is correct.
        // This test mainly focuses on if the evaluation of h1 are correct.
        let mut rng = ark_std::test_rng();

        let hyrax_params = HyraxGenerators::<Projective>::setup(&mut rng, 4);
        let poseidon_config = poseidon_canonical_config::<Fr>();
        let mut transcript_p = PoseidonSponge::<Fr>::new(&poseidon_config);
        let mut transcript_v = PoseidonSponge::<Fr>::new(&poseidon_config);

        let three = Fr::one() + Fr::one() + Fr::one();
        let four = three + Fr::one();
        let five = four + Fr::one();
        let six = five + Fr::one();
        let W_i = MatrixWitness {
            A: Matrix::sparse_from_vec(vec![(0, three), (3, six)], 2, 2).unwrap(),
            B: Matrix::sparse_from_vec(vec![(1, four), (2, five)], 2, 2).unwrap(),
            C: Matrix::sparse_from_vec(vec![(0, three), (3, six)], 2, 2).unwrap(),
            E: Matrix::sparse_from_vec(vec![(1, four), (2, five)], 2, 2).unwrap(),
        };
        let rE = (0..log2(W_i.E.len())).map(|_| Fr::rand(&mut rng)).collect();
        let U_i = MatrixWitness::commit(&W_i, &hyrax_params, rE)?;

        let w_i = MatrixWitness {
            A: Matrix::sparse_from_vec(vec![(0, three), (3, six)], 2, 2).unwrap(),
            B: Matrix::sparse_from_vec(vec![(1, four), (2, five)], 2, 2).unwrap(),
            C: Matrix::sparse_from_vec(vec![(0, three), (3, six)], 2, 2).unwrap(),
            E: Matrix::sparse_from_vec(vec![(1, four), (2, five)], 2, 2).unwrap(),
        };
        let rE = (0..log2(W_i.E.len())).map(|_| Fr::rand(&mut rng)).collect();
        let u_i = MatrixWitness::commit(&w_i, &hyrax_params, rE)?;

        let (proof, claim) = PointVsLineMatrix::prove(&mut transcript_p, None, &u_i, &W_i, &w_i)?;

        let result = PointVsLineMatrix::verify(
            &mut transcript_v,
            None,
            &u_i,
            &proof,
            None,
            &claim.mleE2_prime,
            &claim.rE_prime,
        );

        assert!(result.is_ok(), "Verification failed");
        // Check if the re_prime is the same
        let re_verified = result.unwrap();
        assert!(re_verified == claim.rE_prime);
        let mut transcript_v = PoseidonSponge::<Fr>::new(&poseidon_config);

        // Pass the wrong committed instance which should result in a wrong evaluation in h returning an error
        let result = PointVsLineMatrix::verify(
            &mut transcript_v,
            None,
            &U_i,
            &proof,
            None,
            &claim.mleE2_prime,
            &claim.rE_prime,
        );

        assert!(result.is_err(), "Verification was okay when it should fail");

        Ok(())
    }

    #[test]
    fn test_compute_h2_compare() {
        use ark_std::test_rng;

        // Test Case 1: Simple case with sparse pattern
        {
            let vanilla_dense = DenseMultilinearExtension::from_evaluations_slice(
                3,
                &[
                    Fr::zero(),
                    Fr::zero(),
                    Fr::one(),
                    Fr::one(),
                    Fr::zero(),
                    Fr::zero(),
                    Fr::zero(),
                    Fr::one(),
                ],
            );
            let mle_dense = MultilinearExtension::DenseMLE(vanilla_dense.clone());
            let mle_sparse =
                MultilinearExtension::SparseMLE(SparseMultilinearExtension::from_evaluations(
                    3,
                    &[(2, Fr::one()), (3, Fr::one()), (7, Fr::one())],
                ));

            let r0 = [Fr::from(1), Fr::from(2), Fr::from(3)];
            let r1 = [Fr::from(5), Fr::from(6), Fr::from(7)];
            let r1_sub_r0: Vec<Fr> = r1.iter().zip(&r0).map(|(&x, y)| x - y).collect();

            let result_h = compute_h(&vanilla_dense, &r0, &r1_sub_r0).unwrap();
            let result_h2_dense = compute_h2(&mle_dense, &r0, &r1_sub_r0).unwrap();
            let result_h2_sparse = compute_h2(&mle_sparse, &r0, &r1_sub_r0).unwrap();

            assert_eq!(result_h2_dense, result_h2_sparse);
            assert_eq!(result_h2_dense.coeffs(), result_h.coeffs());
            assert_eq!(result_h2_sparse.coeffs(), result_h.coeffs());
        }

        // Test Case 2: Larger size with random values (4 variables)
        {
            let mut rng = test_rng();
            let evaluations: Vec<Fr> = (0..16).map(|_| Fr::rand(&mut rng)).collect();

            let vanilla_dense =
                DenseMultilinearExtension::from_evaluations_vec(4, evaluations.clone());
            let mle_dense = MultilinearExtension::DenseMLE(vanilla_dense.clone());

            // Create sparse version by filtering out small values
            let sparse_evals: Vec<(usize, Fr)> = evaluations
                .iter()
                .enumerate()
                .filter(|(_, &val)| !val.is_zero())
                .map(|(i, &val)| (i, val))
                .collect();
            let mle_sparse = MultilinearExtension::SparseMLE(
                SparseMultilinearExtension::from_evaluations(4, &sparse_evals),
            );

            let r0: Vec<Fr> = (0..4).map(|_| Fr::rand(&mut rng)).collect();
            let r1: Vec<Fr> = (0..4).map(|_| Fr::rand(&mut rng)).collect();
            let r1_sub_r0: Vec<Fr> = r1.iter().zip(&r0).map(|(&x, y)| x - y).collect();

            let result_h = compute_h(&vanilla_dense, &r0, &r1_sub_r0).unwrap();
            let result_h2_dense = compute_h2(&mle_dense, &r0, &r1_sub_r0).unwrap();
            let result_h2_sparse = compute_h2(&mle_sparse, &r0, &r1_sub_r0).unwrap();

            assert_eq!(
                result_h2_dense, result_h2_sparse,
                "Random 4-var case: dense vs sparse mismatch"
            );
            assert_eq!(
                result_h2_dense.coeffs(),
                result_h.coeffs(),
                "Random 4-var case: h2_dense vs h mismatch"
            );
        }

        // Test Case 3: Edge case - all zeros except one
        {
            let mut evaluations = vec![Fr::zero(); 8];
            evaluations[5] = Fr::from(42);

            let vanilla_dense = DenseMultilinearExtension::from_evaluations_vec(3, evaluations);
            let mle_dense = MultilinearExtension::DenseMLE(vanilla_dense.clone());
            let mle_sparse = MultilinearExtension::SparseMLE(
                SparseMultilinearExtension::from_evaluations(3, &[(5, Fr::from(42))]),
            );

            let r0 = [Fr::from(7), Fr::from(11), Fr::from(13)];
            let r1 = [Fr::from(17), Fr::from(19), Fr::from(23)];
            let r1_sub_r0: Vec<Fr> = r1.iter().zip(&r0).map(|(&x, y)| x - y).collect();

            let result_h = compute_h(&vanilla_dense, &r0, &r1_sub_r0).unwrap();
            let result_h2_dense = compute_h2(&mle_dense, &r0, &r1_sub_r0).unwrap();
            let result_h2_sparse = compute_h2(&mle_sparse, &r0, &r1_sub_r0).unwrap();

            assert_eq!(
                result_h2_dense, result_h2_sparse,
                "Single non-zero case: dense vs sparse mismatch"
            );
            assert_eq!(
                result_h2_dense.coeffs(),
                result_h.coeffs(),
                "Single non-zero case: h2_dense vs h mismatch"
            );
        }

        // Test Case 4: Edge case - all ones (dense case)
        {
            let evaluations = vec![Fr::one(); 16];

            let vanilla_dense =
                DenseMultilinearExtension::from_evaluations_vec(4, evaluations.clone());
            let mle_dense = MultilinearExtension::DenseMLE(vanilla_dense.clone());
            let sparse_evals: Vec<(usize, Fr)> = (0..16).map(|i| (i, Fr::one())).collect();
            let mle_sparse = MultilinearExtension::SparseMLE(
                SparseMultilinearExtension::from_evaluations(4, &sparse_evals),
            );

            let r0 = [Fr::from(2), Fr::from(3), Fr::from(5), Fr::from(7)];
            let r1 = [Fr::from(11), Fr::from(13), Fr::from(17), Fr::from(19)];
            let r1_sub_r0: Vec<Fr> = r1.iter().zip(&r0).map(|(&x, y)| x - y).collect();

            let result_h = compute_h(&vanilla_dense, &r0, &r1_sub_r0).unwrap();
            let result_h2_dense = compute_h2(&mle_dense, &r0, &r1_sub_r0).unwrap();
            let result_h2_sparse = compute_h2(&mle_sparse, &r0, &r1_sub_r0).unwrap();

            assert_eq!(
                result_h2_dense, result_h2_sparse,
                "All ones case: dense vs sparse mismatch"
            );
            assert_eq!(
                result_h2_dense.coeffs(),
                result_h.coeffs(),
                "All ones case: h2_dense vs h mismatch"
            );
        }

        // Test Case 5: Alternating pattern
        {
            let evaluations: Vec<Fr> = (0..32)
                .map(|i| {
                    if i % 2 == 0 {
                        Fr::from(i as u64 + 1)
                    } else {
                        Fr::zero()
                    }
                })
                .collect();

            let vanilla_dense =
                DenseMultilinearExtension::from_evaluations_vec(5, evaluations.clone());
            let mle_dense = MultilinearExtension::DenseMLE(vanilla_dense.clone());

            let sparse_evals: Vec<(usize, Fr)> = evaluations
                .iter()
                .enumerate()
                .filter(|(_, &val)| !val.is_zero())
                .map(|(i, &val)| (i, val))
                .collect();
            let mle_sparse = MultilinearExtension::SparseMLE(
                SparseMultilinearExtension::from_evaluations(5, &sparse_evals),
            );

            let r0 = [
                Fr::from(1),
                Fr::from(4),
                Fr::from(9),
                Fr::from(16),
                Fr::from(25),
            ];
            let r1 = [
                Fr::from(36),
                Fr::from(49),
                Fr::from(64),
                Fr::from(81),
                Fr::from(100),
            ];
            let r1_sub_r0: Vec<Fr> = r1.iter().zip(&r0).map(|(&x, y)| x - y).collect();

            let result_h = compute_h(&vanilla_dense, &r0, &r1_sub_r0).unwrap();
            let result_h2_dense = compute_h2(&mle_dense, &r0, &r1_sub_r0).unwrap();
            let result_h2_sparse = compute_h2(&mle_sparse, &r0, &r1_sub_r0).unwrap();

            assert_eq!(
                result_h2_dense, result_h2_sparse,
                "Alternating pattern case: dense vs sparse mismatch"
            );
            assert_eq!(
                result_h2_dense.coeffs(),
                result_h.coeffs(),
                "Alternating pattern case: h2_dense vs h mismatch"
            );
        }

        // Test Case 6: Very sparse case (only corner evaluations)
        {
            let mut evaluations = vec![Fr::zero(); 16];
            evaluations[0] = Fr::from(10);
            evaluations[15] = Fr::from(20);

            let vanilla_dense = DenseMultilinearExtension::from_evaluations_vec(4, evaluations);
            let mle_dense = MultilinearExtension::DenseMLE(vanilla_dense.clone());
            let mle_sparse =
                MultilinearExtension::SparseMLE(SparseMultilinearExtension::from_evaluations(
                    4,
                    &[(0, Fr::from(10)), (15, Fr::from(20))],
                ));

            let r0 = [Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero()];
            let r1 = [Fr::one(), Fr::one(), Fr::one(), Fr::one()];
            let r1_sub_r0: Vec<Fr> = r1.iter().zip(&r0).map(|(&x, y)| x - y).collect();

            let result_h = compute_h(&vanilla_dense, &r0, &r1_sub_r0).unwrap();
            let result_h2_dense = compute_h2(&mle_dense, &r0, &r1_sub_r0).unwrap();
            let result_h2_sparse = compute_h2(&mle_sparse, &r0, &r1_sub_r0).unwrap();

            assert_eq!(
                result_h2_dense, result_h2_sparse,
                "Corner values case: dense vs sparse mismatch"
            );
            assert_eq!(
                result_h2_dense.coeffs(),
                result_h.coeffs(),
                "Corner values case: h2_dense vs h mismatch"
            );

            // Verify the interpolation property: h(0) should be MLE(r0) and h(1) should be MLE(r1)
            assert_eq!(
                result_h.evaluate(&Fr::zero()),
                Fr::from(10),
                "h(0) should equal MLE(r0)"
            );
            assert_eq!(
                result_h.evaluate(&Fr::one()),
                Fr::from(20),
                "h(1) should equal MLE(r1)"
            );
        }
    }
}
