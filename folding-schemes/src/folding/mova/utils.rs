use std::sync::Arc;
use std::time::Instant;

use ark_ec::CurveGroup;
use ark_ff::PrimeField;
use ark_std::log2;
use ark_std::One;

use crate::utils::mle::dense_vec_to_dense_mle;
use crate::utils::virtual_polynomial::{build_eq_x_r_vec, eq_eval, VirtualPolynomial};
use crate::Error;

use super::{CommittedInstance, Witness};

pub fn compute_c<F: PrimeField>(
    mleE1_prime: F,
    mleE2_prime: F,
    beta: F,
    rE1: &[F],
    rE2: &[F],
    rE_prime: &[F],
) -> Result<F, Error> {
    Ok(mleE1_prime * eq_eval(rE1, rE_prime)? + beta * mleE2_prime * eq_eval(rE2, rE_prime)?)
}

pub fn compute_g<C: CurveGroup>(
    ci1: &CommittedInstance<C>,
    ci2: &CommittedInstance<C>,
    w1: &Witness<C>,
    w2: &Witness<C>,
    beta: &C::ScalarField,
) -> Result<VirtualPolynomial<C::ScalarField>, Error>
where
    C::ScalarField: PrimeField,
{
    let time = Instant::now();
    if w1.E.len() != w2.E.len() {
        return Err(Error::NotEqual);
    }

    let vars = log2(w1.E.len()) as usize;


    let mut g = VirtualPolynomial::<C::ScalarField>::new(vars);
    println!("{:?}", time.elapsed());
    let eq_rE1 = build_eq_x_r_vec(&ci1.rE)?;
    let eq_rE2 = build_eq_x_r_vec(&ci2.rE)?;
    println!("{:?}", time.elapsed());

    let eq_rE1_mle = dense_vec_to_dense_mle(vars, &eq_rE1);
    let eq_rE2_mle = dense_vec_to_dense_mle(vars, &eq_rE2);
    println!("{:?}", time.elapsed());


    let mleE1 = dense_vec_to_dense_mle(log2(w1.E.len()) as usize, &w1.E);
    let mleE2 = dense_vec_to_dense_mle(log2(w2.E.len()) as usize, &w2.E);
    println!("{:?}", time.elapsed());


    g.add_mle_list(
        [Arc::new(mleE1.clone()), Arc::new(eq_rE1_mle.clone())],
        C::ScalarField::one(),
    )?;
    g.add_mle_list(
        [Arc::new(mleE2.clone()), Arc::new(eq_rE2_mle.clone())],
        *beta,
    )?;

    Ok(g)
}
