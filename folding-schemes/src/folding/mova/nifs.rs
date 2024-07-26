use ark_crypto_primitives::sponge::Absorb;
use ark_ec::{CurveGroup, Group};
use ark_ff::PrimeField;

use ark_poly::MultilinearExtension;
use ark_std::{log2, Zero};
use std::marker::PhantomData;
use std::time::Instant;

use super::{CommittedInstance, InstanceWitness, Witness};
use crate::arith::r1cs::R1CS;
use crate::commitment::CommitmentScheme;
use crate::transcript::Transcript;

use crate::utils::mle::dense_vec_to_dense_mle;

use crate::utils::vec::{hadamard, mat_vec_mul, vec_add, vec_scalar_mul, vec_sub};

use crate::Error;
use crate::folding::mova::pointvsline::{PointVsLine, PointvsLineEvaluationClaim, PointVsLineProof};

/// Proof defines a multifolding proof
pub struct Proof<C: CurveGroup, > {
    pub hg_proof: PointVsLineProof<C>,
    pub mleE1_prime: C::ScalarField,
    pub mleE2_prime: C::ScalarField,
    pub mleT: C::ScalarField,
}

/// Implements the Non-Interactive Folding Scheme described in section 4 of
/// [Nova](https://eprint.iacr.org/2021/370.pdf)
pub struct NIFS<C: CurveGroup, CS: CommitmentScheme<C>, T: Transcript<C::ScalarField>, > {
    _c: PhantomData<C>,
    _cp: PhantomData<CS>,
    _ct: PhantomData<T>,
}

impl<C: CurveGroup, CS: CommitmentScheme<C>, T: Transcript<C::ScalarField>, >
NIFS<C, CS, T>
    where
        <C as Group>::ScalarField: Absorb,
{
    // compute_T: compute cross-terms T
    pub fn compute_T(
        r1cs: &R1CS<C::ScalarField>,
        u1: C::ScalarField,
        u2: C::ScalarField,
        z1: &[C::ScalarField],
        z2: &[C::ScalarField],
    ) -> Result<Vec<C::ScalarField>, Error> {
        let (A, B, C) = (r1cs.A.clone(), r1cs.B.clone(), r1cs.C.clone());

        // this is parallelizable (for the future)
        let Az1 = mat_vec_mul(&A, z1)?;
        let Bz1 = mat_vec_mul(&B, z1)?;
        let Cz1 = mat_vec_mul(&C, z1)?;
        let Az2 = mat_vec_mul(&A, z2)?;
        let Bz2 = mat_vec_mul(&B, z2)?;
        let Cz2 = mat_vec_mul(&C, z2)?;

        let Az1_Bz2 = hadamard(&Az1, &Bz2)?;
        let Az2_Bz1 = hadamard(&Az2, &Bz1)?;
        let u1Cz2 = vec_scalar_mul(&Cz2, &u1);
        let u2Cz1 = vec_scalar_mul(&Cz1, &u2);

        vec_sub(&vec_sub(&vec_add(&Az1_Bz2, &Az2_Bz1)?, &u1Cz2)?, &u2Cz1)
    }

    pub fn fold_witness(
        r: C::ScalarField,
        w1: &Witness<C>,
        w2: &Witness<C>,
        T: &[C::ScalarField],
    ) -> Result<Witness<C>, Error> {
        let r2 = r * r;
        let E: Vec<C::ScalarField> = vec_add(
            &vec_add(&w1.E, &vec_scalar_mul(T, &r))?,
            &vec_scalar_mul(&w2.E, &r2),
        )?;
        let W: Vec<C::ScalarField> = w1.W.iter().zip(&w2.W).map(|(a, b)| *a + (r * b)).collect();

        let rW = w1.rW + r * w2.rW;
        Ok(Witness::<C> { E, W, rW })
    }

    pub fn fold_committed_instance(
        rho: C::ScalarField,
        ci1: &CommittedInstance<C>, // U_i
        ci2: &CommittedInstance<C>, // u_i
        rE_prime: &[C::ScalarField],
        mleE1_prime: &C::ScalarField,
        mleE2_prime: &C::ScalarField,
        mleT: &C::ScalarField,
    ) -> Result<CommittedInstance<C>, Error> {
        let r2 = rho * rho;
        let mleE = *mleE1_prime + rho * mleT + r2 * mleE2_prime;
        let u = ci1.u + rho * ci2.u;
        let cmW = ci1.cmW + ci2.cmW.mul(rho);
        let x = ci1
            .x
            .iter()
            .zip(&ci2.x)
            .map(|(a, b)| *a + (rho * b))
            .collect::<Vec<C::ScalarField>>();

        Ok(CommittedInstance::<C> {
            rE: Vec::from(rE_prime),
            mleE,
            u,
            cmW,
            x,
        })
    }

    #[allow(clippy::type_complexity)]
    pub fn prove(
        _cs_prover_params: &CS::ProverParams,
        r1cs: &R1CS<C::ScalarField>,
        transcript: &mut impl Transcript<C::ScalarField>,
        ci1: &CommittedInstance<C>,
        ci2: &CommittedInstance<C>,
        w1: &Witness<C>,
        w2: &Witness<C>,
    ) -> Result<(Proof<C>, InstanceWitness<C>), Error> {
        let start = Instant::now();
        let elapsed = start.elapsed();
        println!("Time before homogenization point-vs-line {:?}", elapsed);
        let (
            hg_proof,
            PointvsLineEvaluationClaim {
                mleE1_prime,
                mleE2_prime,
                rE_prime,
            },
        ) = PointVsLine::<C, T>::prove(transcript, ci1, ci2, w1, w2)?;
        let elapsed = start.elapsed();
        println!("Time after homogenization point-vs-line {:?}", elapsed);


        transcript.absorb(&mleE1_prime);
        transcript.absorb(&mleE2_prime);

        let z1: Vec<C::ScalarField> = [vec![ci1.u], ci1.x.to_vec(), w1.W.to_vec()].concat();
        let z2: Vec<C::ScalarField> = [vec![ci2.u], ci2.x.to_vec(), w2.W.to_vec()].concat();

        let elapsed = start.elapsed();
        println!("Time before computing T {:?}", elapsed);

        let T = Self::compute_T(r1cs, ci1.u, ci2.u, &z1, &z2)?;

        let elapsed = start.elapsed();
        println!("Time after computing T {:?}", elapsed);


        let vars = log2(w1.E.len()) as usize;

        if log2(T.len()) as usize != vars {
            return Err(Error::NotEqual);
        }

        let elapsed = start.elapsed();
        println!("Time before mleT evaluation {:?}", elapsed);
        let mleT = dense_vec_to_dense_mle(vars, &T);
        let mleT_evaluated = mleT.evaluate(&rE_prime).ok_or(Error::EvaluationFail)?;

        let elapsed = start.elapsed();
        println!("Time after mleT evaluation {:?}", elapsed);

        transcript.absorb(&mleT_evaluated);

        let rho_scalar = C::ScalarField::from_le_bytes_mod_order(b"rho");
        transcript.absorb(&rho_scalar);
        let rho: C::ScalarField = transcript.get_challenge();
        let _r2 = rho * rho;

        let elapsed = start.elapsed();
        println!("Time before start folding {:?}", elapsed);
        let fold = Ok((
            Proof::<C> {
                hg_proof,
                mleE1_prime,
                mleE2_prime,
                mleT: mleT_evaluated,
            },
            InstanceWitness {
                ci: Self::fold_committed_instance(
                    rho,
                    ci1,
                    ci2,
                    &rE_prime,
                    &mleE1_prime,
                    &mleE2_prime,
                    &mleT_evaluated,
                )?,
                w: Self::fold_witness(rho, w1, w2, &T)?,
            },
        ));
        let elapsed = start.elapsed();
        println!("Time after folding {:?}", elapsed);
        fold
    }

    /// verify implements NIFS.V logic described in [Nova](https://eprint.iacr.org/2021/370.pdf)'s
    /// section 4. It returns the folded Committed Instance
    pub fn verify(
        transcript: &mut impl Transcript<C::ScalarField>,
        ci1: &CommittedInstance<C>,
        ci2: &CommittedInstance<C>,
        proof: &Proof<C,>
    ) -> Result<CommittedInstance<C>, Error> {
        let rE_prime = PointVsLine::<C, T>::verify(
            transcript,
            ci1,
            ci2,
            &proof.hg_proof,
            &proof.mleE1_prime,
            &proof.mleE2_prime,
        )?;

        transcript.absorb(&proof.mleE1_prime);
        transcript.absorb(&proof.mleE2_prime);
        transcript.absorb(&proof.mleT);

        let rho_scalar = C::ScalarField::from_le_bytes_mod_order(b"rho");
        transcript.absorb(&rho_scalar);
        let rho: C::ScalarField = transcript.get_challenge();

        NIFS::<C, CS, T>::fold_committed_instance(
            rho,
            ci1,
            ci2,
            &rE_prime,
            &proof.mleE1_prime,
            &proof.mleE2_prime,
            &proof.mleT,
        )
    }

    pub fn prove_expansion(
        transcript: &mut impl Transcript<C::ScalarField>,
        ci: &CommittedInstance<C>,
        w: &Witness<C>,
    ) -> Result<CommittedInstance<C>, Error> {
        // accept only R1CS witness
        if w.E.iter().any(|x| !x.is_zero()) {
            return Err(Error::NotEqual);
        }

        let vars = log2(w.E.len()) as usize;

        let rE_scalar = C::ScalarField::from_le_bytes_mod_order(b"rE");
        transcript.absorb(&rE_scalar);
        let rE: Vec<C::ScalarField> = transcript.get_challenges(vars);

        Ok(CommittedInstance {
            rE,
            mleE: C::ScalarField::zero(),
            u: ci.u,
            cmW: ci.cmW,
            x: ci.x.clone(),
        })
    }

    // Do the verifier's expansion part.
    pub fn verify_expansion(
        transcript: &mut impl Transcript<C::ScalarField>,
        ci: &CommittedInstance<C>,
        n_vars_mleE: usize,
    ) -> Result<CommittedInstance<C>, Error> {
        let rE_scalar = C::ScalarField::from_le_bytes_mod_order(b"rE");
        transcript.absorb(&rE_scalar);
        let rE: Vec<C::ScalarField> = transcript.get_challenges(n_vars_mleE);

        Ok(CommittedInstance {
            rE,
            mleE: C::ScalarField::zero(),
            u: ci.u,
            cmW: ci.cmW,
            x: ci.x.clone(),
        })
    }
}
