// This module contains utility functions for common Arkworks curves
// to make them easier to use with the Curve trait.

use ark_ec::short_weierstrass::Projective;

#[cfg(feature = "pallas")]
pub mod pallas {
    use ark_pallas::Projective as PallasProjective;
    
    pub fn pallas_curve() -> PallasProjective {
        PallasProjective::zero()
    }
}

#[cfg(feature = "vesta")]
pub mod vesta {
    use ark_vesta::Projective as VestaProjective;
    
    pub fn vesta_curve() -> VestaProjective {
        VestaProjective::zero()
    }
}

#[cfg(feature = "bn254")]
pub mod bn254 {
    use ark_bn254::G1Projective;
    
    pub fn bn254_curve() -> G1Projective {
        G1Projective::zero()
    }
}