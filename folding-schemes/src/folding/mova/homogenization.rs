use std::fmt::Debug;
use std::marker::PhantomData;

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

pub trait Homogenization<C: CurveGroup, T: Transcript<C>> {
    type Proof: Clone + Debug;

    fn prove(
        transcript: &mut impl Transcript<C>,
        ci1: &CommittedInstance<C>,
        ci2: &CommittedInstance<C>,
        w1: &Witness<C>,
        w2: &Witness<C>,
    ) -> Result<(Self::Proof, HomogeneousEvaluationClaim<C>), Error>;

    fn verify(
        transcript: &mut impl Transcript<C>,
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

#[derive(Clone, Debug, Default)]
pub struct SumCheckHomogenization<C: CurveGroup, T: Transcript<C>> {
    _phantom_C: std::marker::PhantomData<C>,
    _phantom_T: std::marker::PhantomData<T>,
}

#[derive(Clone, Debug)]
pub struct PointVsLineProof<C: CurveGroup> {
    pub h1: DensePolynomial<C::ScalarField>,
    pub h2: DensePolynomial<C::ScalarField>,
}

#[derive(Clone, Debug, Default)]
pub struct PointVsLineHomogenization<C: CurveGroup, T: Transcript<C>> {
    _phantom_C: std::marker::PhantomData<C>,
    _phantom_T: std::marker::PhantomData<T>,
}

impl<C: CurveGroup, T: Transcript<C>> Homogenization<C, T> for SumCheckHomogenization<C, T>
where
    <C as Group>::ScalarField: Absorb,
{
    type Proof = SumCheckProof<C::ScalarField>;

    fn prove(
        transcript: &mut impl Transcript<C>,
        ci1: &CommittedInstance<C>,
        ci2: &CommittedInstance<C>,
        w1: &Witness<C>,
        w2: &Witness<C>,
    ) -> Result<(Self::Proof, HomogeneousEvaluationClaim<C>), Error> {
        let vars = log2(w1.E.len()) as usize;

        let beta_scalar = C::ScalarField::from_le_bytes_mod_order(b"beta");
        transcript.absorb(&beta_scalar);
        let beta: C::ScalarField = transcript.get_challenge();

        let g = compute_g(ci1, ci2, w1, w2, &beta)?;

        let sumcheck_proof = IOPSumCheck::<C, T>::prove(&g, transcript)
            .map_err(|err| Error::SumCheckProveError(err.to_string()))?;

        let rE_prime = sumcheck_proof.point.clone();

        let mleE1 = dense_vec_to_dense_mle(vars, &w1.E);
        let mleE2 = dense_vec_to_dense_mle(vars, &w2.E);

        let mleE1_prime = mleE1.evaluate(&rE_prime).ok_or(Error::EvaluationFail)?;
        let mleE2_prime = mleE2.evaluate(&rE_prime).ok_or(Error::EvaluationFail)?;

        Ok((
            sumcheck_proof,
            HomogeneousEvaluationClaim {
                mleE1_prime,
                mleE2_prime,
                rE_prime,
            },
        ))
    }

    fn verify(
        transcript: &mut impl Transcript<C>,
        ci1: &CommittedInstance<C>,
        ci2: &CommittedInstance<C>,
        proof: &Self::Proof,
        mleE1_prime: &C::ScalarField,
        mleE2_prime: &C::ScalarField,
    ) -> Result<
        Vec<<C>::ScalarField>, // rE=rE1'=rE2'
        Error,
    > {
        let beta_scalar = C::ScalarField::from_le_bytes_mod_order(b"beta");
        transcript.absorb(&beta_scalar);
        let beta: C::ScalarField = transcript.get_challenge();

        let vp_aux_info = VPAuxInfo::<C::ScalarField> {
            max_degree: 2,
            num_variables: ci1.rE.len(),
            phantom: PhantomData::<C::ScalarField>,
        };

        // Step 3: Start verifying the sumcheck
        // First, compute the expected sumcheck sum: \sum gamma^j v_j
        let mut sum_evaluation_claims = ci1.mleE;

        sum_evaluation_claims += beta * ci2.mleE;

        // Verify the interactive part of the sumcheck
        let sumcheck_subclaim =
            IOPSumCheck::<C, T>::verify(sum_evaluation_claims, proof, &vp_aux_info, transcript)
                .map_err(|err| Error::SumCheckVerifyError(err.to_string()))?;

        let rE_prime = sumcheck_subclaim.point.clone();

        let g = compute_c(
            *mleE1_prime,
            *mleE2_prime,
            beta,
            &ci1.rE,
            &ci2.rE,
            &rE_prime,
        )?;

        if g != sumcheck_subclaim.expected_evaluation {
            return Err(Error::NotEqual);
        }

        let g_on_rxprime_from_sumcheck_last_eval = DensePolynomial::from_coefficients_slice(
            &proof.proofs.last().ok_or(Error::Empty)?.coeffs,
        )
        .evaluate(rE_prime.last().ok_or(Error::Empty)?);
        if g_on_rxprime_from_sumcheck_last_eval != g {
            return Err(Error::NotEqual);
        }
        if g_on_rxprime_from_sumcheck_last_eval != sumcheck_subclaim.expected_evaluation {
            return Err(Error::NotEqual);
        }

        Ok(rE_prime)
    }
}

impl<C: CurveGroup, T: Transcript<C>> Homogenization<C, T> for PointVsLineHomogenization<C, T>
where
    <C as Group>::ScalarField: Absorb,
{
    type Proof = PointVsLineProof<C>;

    fn prove(
        transcript: &mut impl Transcript<C>,
        ci1: &CommittedInstance<C>,
        ci2: &CommittedInstance<C>,
        w1: &Witness<C>,
        w2: &Witness<C>,
    ) -> Result<(Self::Proof, HomogeneousEvaluationClaim<C>), Error> {
        let vars = log2(w1.E.len()) as usize;

        let mleE1 = dense_vec_to_dense_mle(vars, &w1.E);
        let mleE2 = dense_vec_to_dense_mle(vars, &w2.E);

        let h1 = compute_h(&mleE1, &ci1.rE, &ci2.rE)?;
        let h2 = compute_h(&mleE2, &ci1.rE, &ci2.rE)?;

        transcript.absorb_vec(h1.coeffs());
        transcript.absorb_vec(h2.coeffs());

        let beta_scalar = C::ScalarField::from_le_bytes_mod_order(b"beta");
        transcript.absorb(&beta_scalar);
        let beta = transcript.get_challenge();

        let mleE1_prime = h1.evaluate(&beta);
        let mleE2_prime = h2.evaluate(&beta);

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
        transcript: &mut impl Transcript<C>,
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

        if proof.h2.evaluate(&C::ScalarField::one()) != ci1.mleE {
            return Err(Error::NotEqual);
        }

        transcript.absorb_vec(proof.h1.coeffs());
        transcript.absorb_vec(proof.h2.coeffs());

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
    use ark_ff::fields::{Fp64, MontBackend, MontConfig};
    use ark_poly::{DenseMultilinearExtension, DenseUVPolynomial};

    use super::compute_h;

    #[derive(MontConfig)]
    #[modulus = "101"]
    #[generator = "3"]
    pub struct FqConfig;
    pub type Fq = Fp64<MontBackend<FqConfig, 1>>;

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
    }
}
