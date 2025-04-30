// This module contains explicit implementations of the Curve trait
// for common Arkworks curves to make them available to external users.

use crate::Curve;
use ark_ec::short_weierstrass::{Projective, SWCurveConfig};
use ark_r1cs_std::{
    fields::fp::FpVar,
    groups::curves::short_weierstrass::ProjectiveVar,
};

#[cfg(feature = "pallas")]
pub mod pallas {
    use super::*;
    use ark_pallas::{Projective as PallasProjective, Config as PallasConfig};
    
    impl Curve for PallasProjective {
        type Var = ProjectiveVar<PallasConfig, FpVar<ark_pallas::Fq>>;
    }
    
    pub fn pallas_curve() -> PallasProjective {
        PallasProjective::zero()
    }
}

#[cfg(feature = "vesta")]
pub mod vesta {
    use super::*;
    use ark_vesta::{Projective as VestaProjective, Config as VestaConfig};
    
    impl Curve for VestaProjective {
        type Var = ProjectiveVar<VestaConfig, FpVar<ark_vesta::Fq>>;
    }
    
    pub fn vesta_curve() -> VestaProjective {
        VestaProjective::zero()
    }
}

#[cfg(feature = "bn254")]
pub mod bn254 {
    use super::*;
    use ark_bn254::{G1Projective, Parameters as BN254Parameters};
    
    impl Curve for G1Projective {
        type Var = ProjectiveVar<BN254Parameters, FpVar<ark_bn254::Fq>>;
    }
    
    pub fn bn254_curve() -> G1Projective {
        G1Projective::zero()
    }
} 