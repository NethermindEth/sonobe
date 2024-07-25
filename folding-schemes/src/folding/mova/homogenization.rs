use std::fmt::Debug;
use std::marker::PhantomData;
use std::time::Instant;

use ark_crypto_primitives::sponge::Absorb;
use ark_ec::{CurveGroup, Group};
use ark_ff::PrimeField;
use ark_poly::univariate::DensePolynomial;
use ark_poly::DenseMultilinearExtension;
use ark_poly::DenseUVPolynomial;
use ark_poly::MultilinearExtension;
use ark_poly::Polynomial;
use ark_std::log2;
use ark_std::One;
use ark_std::Zero;

use crate::folding::mova::utils::compute_c;
use crate::folding::mova::utils::compute_g;
use crate::folding::mova::CommittedInstance;
use crate::folding::mova::Witness;
use crate::transcript::Transcript;
use crate::utils::mle::dense_vec_to_dense_mle;
use crate::utils::sum_check::structs::IOPProof as SumCheckProof;
use crate::utils::sum_check::IOPSumCheck;
use crate::utils::sum_check::SumCheck;
use crate::utils::virtual_polynomial::VPAuxInfo;
use crate::Error;

pub struct HomogeneousEvaluationClaim<C: CurveGroup> {
    pub mleE1_prime: C::ScalarField,
    pub mleE2_prime: C::ScalarField,
    pub rE_prime: Vec<C::ScalarField>,
}

pub trait Homogenization<C: CurveGroup, T: Transcript<C::ScalarField>> {
    type Proof: Clone + Debug;

    fn prove(
        transcript: &mut impl Transcript<C::ScalarField>,
        ci1: &CommittedInstance<C>,
        ci2: &CommittedInstance<C>,
        w1: &Witness<C>,
        w2: &Witness<C>,
    ) -> Result<(Self::Proof, HomogeneousEvaluationClaim<C>), Error>;

    fn verify(
        transcript: &mut impl Transcript<C::ScalarField>,
        ci1: &CommittedInstance<C>,
        ci2: &CommittedInstance<C>,
        proof: &Self::Proof,
        mleE1_prime: &C::ScalarField,
        mleE2_prime: &C::ScalarField,
    ) -> Result<
        Vec<C::ScalarField>, // rE=rE1'=rE2'.
        Error,
    >;
}

#[derive(Clone, Debug)]
pub struct PointVsLineProof<C: CurveGroup> {
    pub h1: DensePolynomial<C::ScalarField>,
    pub h2: DensePolynomial<C::ScalarField>,
}

#[derive(Clone, Debug, Default)]
pub struct PointVsLineHomogenization<C: CurveGroup, T: Transcript<C::ScalarField>> {
    _phantom_C: std::marker::PhantomData<C>,
    _phantom_T: std::marker::PhantomData<T>,
}

impl<C: CurveGroup, T: Transcript<C::ScalarField>> Homogenization<C, T> for PointVsLineHomogenization<C, T>
where
    <C as Group>::ScalarField: Absorb,
{
    type Proof = PointVsLineProof<C>;

    fn prove(
        transcript: &mut impl Transcript<C::ScalarField>,
        ci1: &CommittedInstance<C>,
        ci2: &CommittedInstance<C>,
        w1: &Witness<C>,
        w2: &Witness<C>,
    ) -> Result<(Self::Proof, HomogeneousEvaluationClaim<C>), Error> {
        let vars = log2(w1.E.len()) as usize;

        let mleE1 = dense_vec_to_dense_mle(vars, &w1.E);
        let mleE2 = dense_vec_to_dense_mle(vars, &w2.E);
        let start = Instant::now();
        let elapsed = start.elapsed();
        println!("Time before computing h {:?}", elapsed);

        let h1 = compute_h(&mleE1, &ci1.rE, &ci2.rE)?;
        let h2 = compute_h(&mleE2, &ci1.rE, &ci2.rE)?;

        let elapsed = start.elapsed();
        println!("Time after computing h1 h2 {:?}", elapsed);

        transcript.absorb(&h1.coeffs());
        transcript.absorb(&h2.coeffs());

        let beta_scalar = C::ScalarField::from_le_bytes_mod_order(b"beta");
        transcript.absorb(&beta_scalar);
        let beta = transcript.get_challenge();

        let elapsed = start.elapsed();
        println!("Time before evaluating h {:?}", elapsed);
        let mleE1_prime = h1.evaluate(&beta);
        let mleE2_prime = h2.evaluate(&beta);
        let elapsed = start.elapsed();
        println!("Time after evaluating h {:?}", elapsed);

        let rE_prime = compute_l(&ci1.rE, &ci2.rE, beta)?;

        Ok((
            Self::Proof { h1, h2 },
            HomogeneousEvaluationClaim {
                mleE1_prime,
                mleE2_prime,
                rE_prime,
            },
        ))
    }

    fn verify(
        transcript: &mut impl Transcript<C::ScalarField>,
        ci1: &CommittedInstance<C>,
        ci2: &CommittedInstance<C>,
        proof: &Self::Proof,
        mleE1_prime: &<C>::ScalarField,
        mleE2_prime: &<C>::ScalarField,
    ) -> Result<
        Vec<<C>::ScalarField>, // rE=rE1'=rE2'.
        Error,
    > {
        if proof.h1.evaluate(&C::ScalarField::zero()) != ci1.mleE {
            return Err(Error::NotEqual);
        }

        if proof.h2.evaluate(&C::ScalarField::one()) != ci2.mleE {
            return Err(Error::NotEqual);
        }

        transcript.absorb(&proof.h1.coeffs());
        transcript.absorb(&proof.h2.coeffs());

        let beta_scalar = C::ScalarField::from_le_bytes_mod_order(b"beta");
        transcript.absorb(&beta_scalar);
        let beta = transcript.get_challenge();

        if *mleE1_prime != proof.h1.evaluate(&beta) {
            return Err(Error::NotEqual);
        }

        if *mleE2_prime != proof.h2.evaluate(&beta) {
            return Err(Error::NotEqual);
        }

        let rE_prime = compute_l(&ci1.rE, &ci2.rE, beta)?;

        Ok(rE_prime)
    }
}

// TODO: Test this.
fn compute_l<F: PrimeField>(r0: &[F], r1: &[F], x: F) -> Result<Vec<F>, Error> {
    if r0.len() != r1.len() {
        return Err(Error::NotEqual);
    }

    Ok(r0.iter().zip(r1).map(|(&r0, &r1)| r0 + x * r1).collect())
}

// TODO: This requires thorough testing.
fn compute_h<F: PrimeField>(
    mle: &DenseMultilinearExtension<F>,
    r0: &[F],
    r1: &[F],
) -> Result<DensePolynomial<F>, Error> {
    let nv = mle.num_vars;
    if r0.len() != r1.len() || r0.len() != nv {
        return Err(Error::NotEqual);
    }

    let mut poly: Vec<DensePolynomial<F>> = mle
        .evaluations
        .iter()
        .map(|&x| DensePolynomial::from_coefficients_slice(&[x]))
        .collect();

    // evaluate single variable of partial point from left to right
    for i in 1..nv + 1 {
        let r = DensePolynomial::<F>::from_coefficients_slice(&[r0[i - 1], r1[i - 1] - r0[i - 1]]);
        for b in 0..1 << (nv - i) {
            let left = &poly[b << 1];
            let right = &poly[(b << 1) + 1];
            poly[b] = left + &(&r * &(right - left));
        }
    }

    // Is it even fine?
    Ok(poly.swap_remove(0))
}

#[cfg(test)]
mod tests {
    use super::compute_h;
    use ark_pallas::Fq;
    use ark_poly::{DenseMultilinearExtension, DenseUVPolynomial};

    #[test]
    fn test_compute_h() {
        let mle = DenseMultilinearExtension::from_evaluations_slice(1, &[Fq::from(1), Fq::from(2)]);
        let r0 = [Fq::from(5)];
        let r1 = [Fq::from(6)];
        let result = compute_h(&mle, &r0, &r1).unwrap();
        assert_eq!(
            result,
            DenseUVPolynomial::from_coefficients_slice(&[Fq::from(6), Fq::from(1)])
        );

        let mle = DenseMultilinearExtension::from_evaluations_slice(1, &[Fq::from(1), Fq::from(2)]);
        let r0 = [Fq::from(4)];
        let r1 = [Fq::from(7)];
        let result = compute_h(&mle, &r0, &r1).unwrap();
        assert_eq!(
            result,
            DenseUVPolynomial::from_coefficients_slice(&[Fq::from(5), Fq::from(3)])
        );

        let mle = DenseMultilinearExtension::from_evaluations_slice(
            2,
            &[Fq::from(1), Fq::from(2), Fq::from(3), Fq::from(4)],
        );
        let r0 = [Fq::from(5), Fq::from(4)];
        let r1 = [Fq::from(2), Fq::from(7)];
        let result = compute_h(&mle, &r0, &r1).unwrap();
        assert_eq!(
            result,
            DenseUVPolynomial::from_coefficients_slice(&[Fq::from(14), Fq::from(3)])
        );
        let mle = DenseMultilinearExtension::from_evaluations_slice(
            3,
            &[
                Fq::from(1),
                Fq::from(2),
                Fq::from(3),
                Fq::from(4),
                Fq::from(5),
                Fq::from(6),
                Fq::from(7),
                Fq::from(8),
            ],
        );
        let r0 = [Fq::from(1), Fq::from(2), Fq::from(3)];
        let r1 = [Fq::from(5), Fq::from(6), Fq::from(7)];
        let result = compute_h(&mle, &r0, &r1).unwrap();
        assert_eq!(
            result,
            DenseUVPolynomial::from_coefficients_slice(&[Fq::from(18), Fq::from(28)])
        );
    }
    #[test]
    fn test_compute_h_errors() {
        let mle = DenseMultilinearExtension::from_evaluations_slice(1, &[Fq::from(1), Fq::from(2)]);
        let r0 = [Fq::from(5)];
        let r1 = [];
        let result = compute_h(&mle, &r0, &r1);
        assert!(result.is_err());

        let mle = DenseMultilinearExtension::from_evaluations_slice(
            2,
            &[Fq::from(1), Fq::from(2), Fq::from(1), Fq::from(2)],
        );
        let r0 = [Fq::from(4)];
        let r1 = [Fq::from(7)];
        let result = compute_h(&mle, &r0, &r1);
        assert!(result.is_err())
    }
}
