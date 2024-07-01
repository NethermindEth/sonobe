use std::fmt::Debug;
use std::marker::PhantomData;

use ark_crypto_primitives::sponge::Absorb;
use ark_ec::{CurveGroup, Group};
use ark_ff::PrimeField;
use ark_poly::univariate::DensePolynomial;
use ark_poly::DenseUVPolynomial;
use ark_poly::MultilinearExtension;
use ark_poly::Polynomial;
use ark_std::log2;

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

pub trait Homogenization<C: CurveGroup, T: Transcript<C>> {
    type Proof: Clone + Debug;

    fn prove(
        transcript: &mut impl Transcript<C>,
        ci1: &CommittedInstance<C>,
        ci2: &CommittedInstance<C>,
        w1: &Witness<C>,
        w2: &Witness<C>,
    ) -> Result<
        (
            Self::Proof,
            C::ScalarField,      // mleE1',
            C::ScalarField,      // mleE2',
            Vec<C::ScalarField>, // rE=rE1'=rE2'.
        ),
        Error,
    >;

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
    ) -> Result<
        (
            Self::Proof,
            <C>::ScalarField,      // mleE1',
            <C>::ScalarField,      // mleE2',
            Vec<<C>::ScalarField>, // rE=rE1'=rE2'.
        ),
        Error,
    > {
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

        Ok((sumcheck_proof, mleE1_prime, mleE2_prime, rE_prime))
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
            max_degree: 1,
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
