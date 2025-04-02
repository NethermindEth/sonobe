#![allow(dead_code)]
use crate::commitment::hyrax::{Hyrax, HyraxGenerators};
/// Mova-like folding for matrix multiplications as described in "Folding and Lookup Arguments for Proving Inference of Deep Learning Models" by Nethermind Research
/// Currently, we are not interested in the hiding properties, so we ignore the hiding factors and focus on the succinctness property.
/// Please note the code could be easily extended so support hiding.
use crate::folding::nova::nifs::pointvsline::{
    PointVsLine, PointVsLineEvaluationClaimMatrix, PointVsLineMatrix, PointVsLineProofMatrix,
};
use crate::transcript::Transcript;
use crate::utils::mle::MultilinearExtension;
use crate::utils::vec::is_zero_vec;
use crate::{Curve, Error};
use ark_crypto_primitives::sponge::Absorb;
use ark_ff::PrimeField;
use ark_std::{log2, marker::PhantomData, rand::RngCore, One, UniformRand, Zero};
use matrex::Matrix;
use num_integer::Roots;

/// Represents a relaxed committed relation for matrix multiplication folded instances.
/// Where A, B, C are nxn matrices such that:
/// A * B = u * C + E
/// where u is a scalar, and E is the error matrix.
/// When u=1 and E is the zero matrix, we have the simple committed relation in which A * B = C.
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct RelaxedCommittedRelation<C: Curve> {
    pub cmA: Vec<C>,             // Commitment to matrix A. cmA = commitment(A).
    pub cmB: Vec<C>,             // Commitment to matrix B. cmB = commitment(B).
    pub cmC: Vec<C>,             // Commitment to matrix C. cmC = commitment(C).
    pub u: C::ScalarField,       // Scalar used in the folding.
    pub mleE: C::ScalarField, // v = mle[E](rE) in MOVA notation. Multilinear extension of matrix E evaluated at random point rE.
    pub rE: Vec<C::ScalarField>, // Random point where MLE is evaluated. Size of 2 * log2(n).
}

impl<C: Curve> Absorb for RelaxedCommittedRelation<C> {
    fn to_sponge_bytes(&self, dest: &mut Vec<u8>) {
        C::ScalarField::batch_to_sponge_bytes(&self.to_sponge_field_elements_as_vec(), dest);
    }

    fn to_sponge_field_elements<F: PrimeField>(&self, dest: &mut Vec<F>) {
        for i in 0..self.cmA.len() {
            self.cmA[i].to_native_sponge_field_elements(dest);
            self.cmB[i].to_native_sponge_field_elements(dest);
            self.cmC[i].to_native_sponge_field_elements(dest);
        }
        self.u.to_sponge_field_elements(dest);
        self.mleE.to_sponge_field_elements(dest);
        self.rE.to_sponge_field_elements(dest);
    }
}

impl<C: Curve> RelaxedCommittedRelation<C> {
    /// Checks if a Relaxed Committed Relation is simple (has not been folded).
    fn is_simple(&self) -> bool {
        self.u == C::ScalarField::from(1) && self.mleE == C::ScalarField::zero()
    }

    /// Checks if a Relaxed Committed Relation is accumulated (has been folded).
    fn is_accumulated(&self) -> bool {
        self.u != C::ScalarField::from(1) && self.mleE != C::ScalarField::zero()
    }
}
/// Represents the private inputs for the protocol (witness)
/// A, B, C, E are matrices such that A * B = u* C + E
#[derive(Debug, Clone, PartialEq)]
pub struct Witness<C: Curve> {
    pub A: Matrix<C::ScalarField>,
    pub B: Matrix<C::ScalarField>,
    pub C: Matrix<C::ScalarField>,
    pub E: Matrix<C::ScalarField>,
}

impl<C: Curve> Witness<C> {
    pub fn new<const H: bool>(
        a: Matrix<C::ScalarField>,
        b: Matrix<C::ScalarField>,
        c: Matrix<C::ScalarField>,
        e: Matrix<C::ScalarField>,
    ) -> Self {
        Self {
            A: a,
            B: b,
            C: c,
            E: e,
        }
    }

    /// Commits to a witness W and produces a RelaxedCommittedRelation
    /// # Parameters
    /// * `self` - Witness instance to be committed.
    /// * `params` - Commitment scheme parameters.
    /// * `rE` - Random evaluation point for the committed instance.
    pub fn commit(
        &self,
        params: &HyraxGenerators<C>,
        rE: Vec<C::ScalarField>,
    ) -> Result<RelaxedCommittedRelation<C>, Error> {
        let mut mleE = C::ScalarField::zero();
        if !Matrix::is_empty(&self.E) {
            let mle = MultilinearExtension::from_evaluations(&self.E, log2(self.E.len()) as usize);
            mleE = mle.evaluate(&rE);
        }
        // Simple Dense Commit
        // let com_a = Hyrax::commit(self.A.as_dense_slice().unwrap(), params)?;
        // let com_b = Hyrax::commit(self.B.as_dense_slice().unwrap(), params)?;
        // let com_c = Hyrax::commit(self.C.as_dense_slice().unwrap(), params)?;

        // Sparse commit
        let com_a = Hyrax::commit_sparse_matrix(self.A.as_sparse_slice().unwrap(), params)?;
        let com_b = Hyrax::commit_sparse_matrix(self.B.as_sparse_slice().unwrap(), params)?;
        let com_c = if self.C.is_dense() {
            Hyrax::commit(self.C.as_dense_slice().unwrap(), params)?
        } else {
            Hyrax::commit_sparse_matrix(self.C.as_sparse_slice().unwrap(), params)?
        };

        // Batch dense commit
        // let a_slice = self.A.as_dense_slice().unwrap();
        // let b_slice = self.B.as_dense_slice().unwrap();
        // let c_slice = self.C.as_dense_slice().unwrap();
        //
        // let slices = [a_slice, b_slice, c_slice];
        //
        // let mut all_commitments = Hyrax::batch_commit(&slices, params)?;
        //
        // let mut commitments_iter = all_commitments.into_iter();
        // let com_a = commitments_iter.next().expect("Expected commitments for A");
        // let com_b = commitments_iter.next().expect("Expected commitments for B");
        // let com_c = commitments_iter.next().expect("Expected commitments for C");

        Ok(RelaxedCommittedRelation {
            cmA: com_a,
            cmB: com_b,
            cmC: com_c,
            u: C::ScalarField::one(),
            mleE,
            rE,
        })
    }
}

/// Represents proof of the matrix folding
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct Proof<C: Curve> {
    /// Proof of the PointVsLine protocol that reduces the number of MLE evaluations.
    pub h_proof: PointVsLineProofMatrix<C>,
    /// Evaluation of the MLE[h2] in the random challenge beta
    pub mleE2_prime: C::ScalarField,
    /// Evaluation of the cross term T in r_acc = MLE[T](r_acc)
    pub mleT: C::ScalarField,
    /// Evaluation of the h polynomial in the random challenge beta
    pub rE_prime: Vec<C::ScalarField>,
}

/// Implements the Non-Interactive Folding Scheme described in section 3.2 of the previous referenced article.
/// `H` specifies whether the NIFS will use a blinding factor
#[allow(clippy::upper_case_acronyms)]
pub struct NIFS<C: Curve, T: Transcript<C::ScalarField>, const H: bool = false> {
    _c: PhantomData<C>,
    _ct: PhantomData<T>,
}

impl<C: Curve, T: Transcript<C::ScalarField>, const H: bool> NIFS<C, T, H> {
    fn new_witness(
        a: Matrix<C::ScalarField>,
        b: Matrix<C::ScalarField>,
        c: Matrix<C::ScalarField>,
        e: Matrix<C::ScalarField>,
    ) -> Witness<C> {
        Witness::new::<H>(a, b, c, e)
    }

    fn new_instance(
        mut rng: impl RngCore,
        params: &HyraxGenerators<C>,
        witness: Witness<C>,
        aux: Vec<C::ScalarField>, // rE in MOVA notation.
    ) -> Result<RelaxedCommittedRelation<C>, Error> {
        let mut rE = aux;
        if is_zero_vec(&rE) {
            // means that we're in a fresh instance, so generate value of length 2 * log2(n)
            rE = (0..2 * log2(witness.E.len().sqrt()))
                .map(|_| C::ScalarField::rand(&mut rng))
                .collect();
        }
        witness.commit(params, rE)
    }

    // Protocol 5 - point 8 (Page 25)
    fn fold_witness(
        alpha: C::ScalarField,           // Random challenge
        simple_wit: &mut Witness<C>,     // Simple witness
        acc_wit: &Witness<C>,            // Accumulated witness
        mut aux: Matrix<C::ScalarField>, // T in Mova's notation
    ) -> Result<Witness<C>, Error> {
        simple_wit.A *= alpha;
        simple_wit.B *= alpha;
        simple_wit.C *= alpha;
        aux *= alpha;

        let a_acc = (&simple_wit.A + &acc_wit.A).unwrap();
        let b_acc = (&simple_wit.B + &acc_wit.B).unwrap();
        let c_acc = (&simple_wit.C + &acc_wit.C).unwrap();
        let e_acc = (aux + &acc_wit.E).unwrap();

        Ok(Witness::<C> {
            A: a_acc,
            B: b_acc,
            C: c_acc,
            E: e_acc,
        })
    }

    /// Protocol 5 for MOVA-like matrix folding
    /// Returns a proof for the pt-vs-line operations along with the folded committed instance
    /// instances and witness
    #[allow(clippy::type_complexity)]
    pub fn prove(
        transcript: &mut T,
        pp_hash: C::ScalarField,
        simple_witness: &mut Witness<C>,
        simple_instance: &RelaxedCommittedRelation<C>,
        acc_witness: &Witness<C>,
        acc_instance: &RelaxedCommittedRelation<C>,
    ) -> Result<(Witness<C>, RelaxedCommittedRelation<C>, Proof<C>), Error> {
        // Verify instances have the correct form.
        // 2 simple instances can be folded, a simple and an accumulated instance can also be folded. 2 accumulated instances cannot be folded
        if simple_instance.is_accumulated() {
            return if acc_instance.is_simple() {
                Err(Error::Other(String::from(
                    "Parameters were passed in the wrong order. They need to be reordered.",
                )))
            } else {
                Err(Error::Other(String::from(
                    "Cannot fold 2 accumulated instances.",
                )))
            };
        }

        transcript.absorb(&pp_hash);
        transcript.absorb(simple_instance);
        transcript.absorb(acc_instance);

        // Protocol 5 - Steps 1-3
        let (
            h_proof,
            PointVsLineEvaluationClaimMatrix {
                mleE2_prime,
                rE_prime,
            },
        ) = PointVsLineMatrix::<C, T>::prove(
            transcript,
            None,
            acc_instance,
            simple_witness,
            acc_witness,
        )?;

        transcript.absorb(&mleE2_prime);

        // Compute cross term T
        let mut A1B2 = (&simple_witness.A * &acc_witness.B).unwrap(); // Representing A1B2A2B1 later

        let B1A2 = (&acc_witness.A * &simple_witness.B).unwrap();
        A1B2 += B1A2;
        let u2c1 = simple_witness.C.clone() * acc_instance.u;
        let T: Matrix<C::ScalarField> = ((A1B2 - &acc_witness.C).unwrap() - u2c1).unwrap();

        // Compute MLE_T
        let n_vars: usize = log2(simple_witness.E.len()) as usize;
        if log2(T.len()) as usize != n_vars {
            return Err(Error::NotExpectedLength(T.len(), n_vars));
        }

        let mle = MultilinearExtension::from_evaluations(&T, n_vars);
        let mleT_evaluated = mle.evaluate(&rE_prime);

        // Derive alpha
        transcript.absorb(&mleT_evaluated);
        let alpha: C::ScalarField = transcript.get_challenge();

        let ci = Self::fold_committed_instance(
            alpha,
            simple_instance,
            acc_instance,
            &rE_prime,
            &mleE2_prime,
            &mleT_evaluated,
        )?;
        let w = Self::fold_witness(alpha, simple_witness, acc_witness, T)?;

        let proof = Proof::<C> {
            h_proof,
            mleE2_prime,
            mleT: mleT_evaluated,
            rE_prime,
        };
        Ok((w, ci, proof))
    }

    /// It verifies the results from the proof
    /// Both the folding and the pt-vs-line proof
    /// returns the folded committed instance
    pub fn verify(
        transcript: &mut T,
        pp_hash: C::ScalarField,
        simple_instance: &RelaxedCommittedRelation<C>,
        acc_instance: &RelaxedCommittedRelation<C>,
        proof: &Proof<C>,
    ) -> Result<(RelaxedCommittedRelation<C>, Vec<bool>), Error> {
        transcript.absorb(&pp_hash);
        transcript.absorb(simple_instance);
        transcript.absorb(acc_instance);

        // Verify rE_prime
        let rE_prime = PointVsLineMatrix::<C, T>::verify(
            transcript,
            None,
            acc_instance,
            &proof.h_proof,
            None,
            &proof.mleE2_prime,
            &proof.rE_prime,
        )?;

        // Derive alpha
        transcript.absorb(&proof.mleE2_prime);
        transcript.absorb(&proof.mleT);
        let alpha: C::ScalarField = transcript.get_challenge();

        Ok((
            Self::fold_committed_instance(
                alpha,
                simple_instance,
                acc_instance,
                &rE_prime,
                &proof.mleE2_prime,
                &proof.mleT,
            )?,
            vec![],
        ))
    }
}

impl<C: Curve, T: Transcript<C::ScalarField>, const H: bool> NIFS<C, T, H> {
    /// Folds two committed instances into a single one using the provided parameters
    /// Protocol 5 - Step 7 (page 25) for folding instances
    /// # Parameters
    /// * `alpha` - Random challenge used for folding
    /// * `simple_instance` - The simple (unfolded) instance
    /// * `acc_instance` - The accumulated (previously folded) instance  
    /// * `rE_prime` - New random evaluation point
    /// * `mleE2_prime` - Evaluation of MLE[E2] at rE_prime
    /// * `mleT` - Evaluation of the crossterm T
    fn fold_committed_instance(
        alpha: C::ScalarField,
        simple_instance: &RelaxedCommittedRelation<C>,
        acc_instance: &RelaxedCommittedRelation<C>,
        rE_prime: &[C::ScalarField],
        mleE2_prime: &C::ScalarField, // v' in Protocol 5
        mleT: &C::ScalarField,
    ) -> Result<RelaxedCommittedRelation<C>, Error> {
        // Step 7
        // Accumulate commitments
        // let com_a_acc = simple_instance.cmA.mul(alpha) + acc_instance.cmA;
        // let com_b_acc = simple_instance.cmB.mul(alpha) + acc_instance.cmB;
        // let com_c_acc = simple_instance.cmC.mul(alpha) + acc_instance.cmC;

        let com_a_acc: Vec<_> = simple_instance
            .cmA
            .iter()
            .zip(&acc_instance.cmA)
            .map(|(a, b)| a.mul(alpha).add(b))
            .collect();

        let com_b_acc: Vec<_> = simple_instance
            .cmB
            .iter()
            .zip(&acc_instance.cmB)
            .map(|(a, b)| a.mul(alpha).add(b))
            .collect();

        let com_c_acc: Vec<_> = simple_instance
            .cmC
            .iter()
            .zip(&acc_instance.cmC)
            .map(|(a, b)| a.mul(alpha).add(b))
            .collect();

        // Update scalars
        let u_acc = alpha + acc_instance.u;
        let mleE = *mleE2_prime + alpha * *mleT;

        Ok(RelaxedCommittedRelation::<C> {
            cmA: com_a_acc,
            cmB: com_b_acc,
            cmC: com_c_acc,
            u: u_acc,
            mleE,
            rE: rE_prime.to_vec(),
        })
    }

    /// Checks the relationship of the accumulated instances and witness
    /// AB == uC + E
    /// commit(M) = storedCommitment of M (for A, B and C)
    /// MLE[E](r) == v
    /// A, B and C are nxn matrices - to be added when we have correct matrices
    fn check_relation(
        witness: &Witness<C>,
        instance: &RelaxedCommittedRelation<C>,
        params: &HyraxGenerators<C>,
    ) -> Result<(), Error> {
        let mut dense_a = witness.A.clone();
        let mut dense_b = witness.B.clone();
        let mut dense_c = witness.C.clone();
        dense_a.to_dense();
        dense_b.to_dense();
        dense_c.to_dense();
        let com_a = Hyrax::commit(dense_a.as_dense_slice().unwrap(), params)?;
        let com_b = Hyrax::commit(dense_b.as_dense_slice().unwrap(), params)?;
        let com_c = Hyrax::commit(dense_c.as_dense_slice().unwrap(), params)?;

        let AB = (&witness.A * &witness.B).unwrap();

        let uc = witness.C.clone() * instance.u;

        let e = (&AB - &uc).unwrap();

        if witness.E != e {
            return Err(Error::NotSatisfied);
        }

        let E = MultilinearExtension::from_evaluations(&witness.E, log2(witness.E.len()) as usize);
        let mleE = E.evaluate(&instance.rE);
        if instance.mleE != mleE {
            return Err(Error::NotSatisfied);
        }

        if instance.cmA != com_a {
            return Err(Error::NotSatisfied);
        }
        if instance.cmB != com_b {
            return Err(Error::NotSatisfied);
        }
        if instance.cmC != com_c {
            return Err(Error::NotSatisfied);
        }
        Ok(())
    }
}

#[cfg(test)]
pub mod tests {
    use super::*;
    use crate::transcript::poseidon::poseidon_canonical_config;
    use ark_crypto_primitives::sponge::{poseidon::PoseidonSponge, CryptographicSponge};
    use ark_pallas::{Fr, Projective};
    use matrex::Matrix;
    use rand::Rng;

    fn random_sparse_matrix<C: Curve>(n: usize, rng: &mut impl RngCore) -> Matrix<C::ScalarField> {
        let elements = (0..n)
            .map(|row| {
                (
                    row * n + rand::thread_rng().gen_range(0..n),
                    C::ScalarField::rand(rng),
                )
            })
            .collect();
        Matrix::sparse_from_vec(elements, n, n).unwrap()
    }

    // Helper functions
    fn get_instances<C: Curve>(
        num: usize,
        n: usize,
        rng: &mut impl RngCore,
        params: &HyraxGenerators<C>,
    ) -> Vec<(Witness<C>, RelaxedCommittedRelation<C>)> {
        (0..num)
            .map(|_| -> (Witness<C>, RelaxedCommittedRelation<C>) {
                // A matrix
                let mut a = random_sparse_matrix::<C>(n, rng);
                // a.to_dense();
                // B matrix
                let mut b = random_sparse_matrix::<C>(n, rng);
                // b.to_dense();
                // C = A * B matrix
                let c = (&a * &b).unwrap();
                // Error matrix initialized to 0s
                let e = Matrix::zero(n, n);

                // Random challenge
                let rE = (0..2 * log2(n))
                    .map(|_| C::ScalarField::rand(rng))
                    .collect();
                // Witness
                let witness = Witness::new::<false>(a, b, c, e);
                let instance = witness.commit(params, rE).unwrap();
                (witness, instance)
            })
            .collect()
    }
    #[test]
    fn test_nifs_mova_matrix_single_fold() {
        // Set up test instances
        let mut rng = ark_std::test_rng();
        let n_instances = 2;
        let mat_dim = 4; // 4x4 matrices

        // Set up transcript and commitment scheme
        let hyrax_params =
            HyraxGenerators::<Projective>::setup(&mut rng, log2(mat_dim * mat_dim) as usize);
        let poseidon_config = poseidon_canonical_config::<Fr>();
        let mut transcript_p = PoseidonSponge::<Fr>::new(&poseidon_config);
        let mut transcript_v = PoseidonSponge::<Fr>::new(&poseidon_config);
        let pp_hash = Fr::rand(&mut rng);

        let mut instances: Vec<(Witness<Projective>, RelaxedCommittedRelation<Projective>)> =
            get_instances::<Projective>(n_instances, mat_dim, &mut rng, &hyrax_params);

        let (left, right) = instances.split_at_mut(1);
        let (ref mut simple_w, ref simple_i) = &mut left[0];
        let (ref acc_w, ref acc_i) = &mut right[0];

        let (wit_acc, inst_acc, proof) = NIFS::<Projective, PoseidonSponge<Fr>>::prove(
            &mut transcript_p,
            pp_hash,
            simple_w,
            &simple_i,
            &acc_w,
            &acc_i,
        )
        .unwrap();

        // Verify
        let (ci_verify, _) = NIFS::<Projective, PoseidonSponge<Fr>>::verify(
            &mut transcript_v,
            pp_hash,
            &simple_i,
            &acc_i,
            &proof,
        )
        .unwrap();

        // Ensure they match
        assert_eq!(inst_acc, ci_verify);

        NIFS::<Projective, PoseidonSponge<Fr>>::check_relation(&wit_acc, &inst_acc, &hyrax_params)
            .expect("Relationship check failed");
    }

    #[test]
    fn test_nifs_mova_matrix_multiple_folds() {
        // Set up test instances
        let mut rng = ark_std::test_rng();
        let n_folds = 10;
        let n_instances = n_folds + 1;
        let mat_dim = 16; // 16x16 matrices

        // Set up transcript and commitment scheme
        let hyrax_params =
            HyraxGenerators::<Projective>::setup(&mut rng, log2(mat_dim * mat_dim) as usize);
        let poseidon_config = poseidon_canonical_config::<Fr>();
        let mut transcript_p = PoseidonSponge::<Fr>::new(&poseidon_config);
        let mut transcript_v = PoseidonSponge::<Fr>::new(&poseidon_config);
        let pp_hash = Fr::rand(&mut rng);

        let mut instances: Vec<(Witness<Projective>, RelaxedCommittedRelation<Projective>)> =
            get_instances::<Projective>(n_instances, mat_dim, &mut rng, &hyrax_params);

        // Keep track of the accumulated state
        let first_instance = instances.remove(0);
        let mut current_acc_wit = first_instance.0;
        let mut current_acc_inst = first_instance.1;

        // Fold through all remaining instances
        for (mut next_w, next_i) in instances {
            // Fold
            let (wit_acc, inst_acc, proof) = NIFS::<Projective, PoseidonSponge<Fr>>::prove(
                &mut transcript_p,
                pp_hash,
                &mut next_w,
                &next_i,
                &current_acc_wit,
                &current_acc_inst,
            )
            .unwrap();

            // Verify
            let (ci_verify, _) = NIFS::<Projective, PoseidonSponge<Fr>>::verify(
                &mut transcript_v,
                pp_hash,
                &next_i,
                &current_acc_inst,
                &proof,
            )
            .unwrap();

            // Ensure they match
            assert_eq!(inst_acc, ci_verify);

            // Update state for next iteration
            NIFS::<Projective, PoseidonSponge<Fr>>::check_relation(
                &wit_acc,
                &inst_acc,
                &hyrax_params,
            )
            .expect("Relationship check failed");

            current_acc_wit = wit_acc;
            current_acc_inst = inst_acc;
        }
    }
}
