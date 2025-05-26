//! Projective implementation of a twisted Hessian curve

use crate::{
    field::Fq,
    ring::RingElement,
};
use core::ops::{
    Add,
    Mul,
    Sub,
};

/// Represents a point [X:Y:Z] in projective coordinates on a twisted Hessian curve
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Projective<const Q: u64> {
    x: RingElement<Q>,
    y: RingElement<Q>,
    z: RingElement<Q>,
}

impl<const Q: u64> Projective<Q> {
    /// Create a new projective point [X:Y:Z]
    pub fn new(x: RingElement<Q>, y: RingElement<Q>, z: RingElement<Q>) -> Self {
        Projective { x, y, z }
    }

    /// Create the identity element [0:-1:1]
    pub fn identity() -> Self {
        let zero = Fq::new(0);
        let one = Fq::new(1);
        let neg_one = Fq::new(Q.checked_sub(1).expect("Q must be gt 1"));

        Projective::new(
            RingElement::from_field(zero),
            RingElement::from_field(neg_one),
            RingElement::from_field(one),
        )
    }

    /// Get the x-coordinate
    pub fn x(&self) -> RingElement<Q> {
        self.x
    }

    /// Get the y-coordinate
    pub fn y(&self) -> RingElement<Q> {
        self.y
    }

    /// Get the z-coordinate
    pub fn z(&self) -> RingElement<Q> {
        self.z
    }

    /// Get the modulus of the underlying field
    pub fn modulus(&self) -> u64 {
        Q
    }

    /// Check if this is the identity point
    pub fn is_identity(&self) -> bool {
        self.x.constant().value() == 0
            && self.y.constant().value() == self.modulus().checked_sub(1).unwrap()
            && self.z.constant().value() == 1
    }

    /// Check if a point is "projectively equal" to another
    pub fn is_equal(&self, other: &Self) -> bool {
        let self_is_zero = self.x.constant().value() == 0
            && self.x.epsilon_coeff().value() == 0
            && self.y.constant().value() == 0
            && self.y.epsilon_coeff().value() == 0
            && self.z.constant().value() == 0
            && self.z.epsilon_coeff().value() == 0;

        let other_is_zero = other.x.constant().value() == 0
            && other.x.epsilon_coeff().value() == 0
            && other.y.constant().value() == 0
            && other.y.epsilon_coeff().value() == 0
            && other.z.constant().value() == 0
            && other.z.epsilon_coeff().value() == 0;

        if self_is_zero || other_is_zero {
            panic!("Attempted equality check with invalid point [0:0:0]");
        }

        // two projective points [X1:Y1:Z1] and [X2:Y2:Z2] are equal if
        // X1*Z2 = X2*Z1 and Y1*Z2 = Y2*Z1 and Z1*X2 = Z2*X1
        let x1z2 = self.x.mul(other.z);
        let x2z1 = other.x.mul(self.z);

        let y1z2 = self.y.mul(other.z);
        let y2z1 = other.y.mul(self.z);

        let z1x2 = self.z.mul(other.x);
        let z2x1 = other.z.mul(self.x);

        x1z2 == x2z1 && y1z2 == y2z1 && z1x2 == z2x1
    }

    /// Check if a point lies on a twisted Hessian curve aX³ + Y³ + Z³ = dXYZ
    pub fn is_on_curve(&self, a: RingElement<Q>, d: RingElement<Q>) -> bool {
        // TODO: maybe we don't need the below check since it's done in curve.rs
        let d_cubed = d.mul(d).mul(d);
        let twenty_seven = RingElement::from_field(Fq::new(27u64.rem_euclid(Q)));
        let twenty_seven_a = twenty_seven.mul(a);
        let term = twenty_seven_a.sub(d_cubed);
        let condition = a.mul(term);

        assert!(
            condition.is_invertible(),
            "Invalid curve parameters: a(27a−d³) must be invertible in the ring"
        );

        // aX³ + Y³ + Z³ = dXYZ
        let x_cubed = self.x.mul(self.x).mul(self.x);
        let y_cubed = self.y.mul(self.y).mul(self.y);
        let z_cubed = self.z.mul(self.z).mul(self.z);

        let axyz = a.mul(x_cubed).add(y_cubed).add(z_cubed);

        let dxyz = d.mul(self.x).mul(self.y).mul(self.z);

        axyz == dxyz
    }

    /// Negate a point: -[X:Y:Z] = [X:Z:Y]
    pub fn negate(&self) -> Self {
        Projective::new(self.x, self.z, self.y)
    }

    /// Add two points on a twisted Hessian curve
    pub fn add(&self, other: &Self, a: RingElement<Q>) -> Self {
        // implementation of Algorithm 3.1 (1) from the paper

        // this is weird though, hessian curve additions are supposed to have a unified formula
        if self.is_equal(other) {
            // X'₃ = Z₂²X₁Z₁ - Y₁²X₂Y₂
            // Y'₃ = Y₂²Y₁Z₁ - aX₁²X₂Z₂
            // Z'₃ = aX₂²X₁Y₁ - Z₁²Y₂Z₂
            let x1_squared = self.x.mul(self.x);
            let y1_squared = self.y.mul(self.y);
            let z1_squared = self.z.mul(self.z);

            let x3 = z1_squared
                .mul(self.x)
                .mul(self.z)
                .sub(y1_squared.mul(self.x).mul(self.y));
            let y3 = y1_squared
                .mul(self.y)
                .mul(self.z)
                .sub(a.mul(x1_squared).mul(self.x).mul(self.z));
            let z3 = a
                .mul(x1_squared)
                .mul(self.x)
                .mul(self.y)
                .sub(z1_squared.mul(self.y).mul(self.z));

            return Projective::new(x3, y3, z3);
        }

        // X₃ = X₁²Y₂Z₂ - X₂²Y₁Z₁
        let x1_squared = self.x.mul(self.x);
        let x2_squared = other.x.mul(other.x);
        let term1 = x1_squared.mul(other.y).mul(other.z);
        let term2 = x2_squared.mul(self.y).mul(self.z);
        let x3 = term1.sub(term2);

        // Y₃ = Z₁²X₂Y₂ - Z₂²X₁Y₁
        let z1_squared = self.z.mul(self.z);
        let z2_squared = other.z.mul(other.z);
        let term3 = z1_squared.mul(other.x).mul(other.y);
        let term4 = z2_squared.mul(self.x).mul(self.y);
        let y3 = term3.sub(term4);

        // Z₃ = Y₁²X₂Z₂ - Y₂²X₁Z₁
        let y1_squared = self.y.mul(self.y);
        let y2_squared = other.y.mul(other.y);
        let term5 = y1_squared.mul(other.x).mul(other.z);
        let term6 = y2_squared.mul(self.x).mul(self.z);
        let z3 = term5.sub(term6);

        // handle special case where both formulas give [0:0:0]
        if x3.constant().value() == 0
            && x3.epsilon_coeff().value() == 0
            && y3.constant().value() == 0
            && y3.epsilon_coeff().value() == 0
            && z3.constant().value() == 0
            && z3.epsilon_coeff().value() == 0
        {
            panic!("Point addition resulted in invalid point [0:0:0]");
        }

        Projective::new(x3, y3, z3)
    }

    /// Double a point on a twisted Hessian curve (specialized point addition)
    pub fn double(&self, a: RingElement<Q>) -> Self {
        self.add(self, a)
    }

    /// Multiply a point by a scalar using double-and-add algorithm
    pub fn scalar_mul(&self, scalar: u64, a: RingElement<Q>) -> Self {
        // TODO: optimize using msm
        let mut result = Projective::identity();
        let mut temp = *self;
        let mut k = scalar;

        while k > 0 {
            if k & 1 == 1 {
                result = result.add(&temp, a);
            }
            temp = temp.double(a);
            k >>= 1;
        }

        result
    }

    /// Verify a & d
    pub fn verify_curve_constraints(a: RingElement<Q>, d: RingElement<Q>) -> bool {
        let twenty_seven = RingElement::from_field(Fq::<Q>::new(27 % Q));
        let twenty_seven_a = twenty_seven.mul(a);

        let d_squared = d.mul(d);
        let d_cubed = d_squared.mul(d);

        let term = twenty_seven_a.sub(d_cubed);
        let condition = a.mul(term);

        condition.is_invertible()
    }
}

#[allow(non_snake_case)]
#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::Fq;

    #[test]
    fn kats_paper_3_1() {
        // using F5[ε] and the curve with a = 1+ε, d = 1+ε as in the paper
        type F5 = Fq<5>;

        // Create field elements
        let field_1 = F5::new(1);
        let field_2 = F5::new(2);
        let field_3 = F5::new(3);
        let field_4 = F5::new(4);

        let a = RingElement::new(field_1, field_1); // 1+ε
        let d = RingElement::new(field_1, field_1); // 1+ε

        let twenty_seven = RingElement::from_field(F5::new(27 % 5)); // 27 mod 5 = 2
        let twenty_seven_a = twenty_seven.mul(a); // 2*(1+ε) = 2+2ε
        let d_squared = d.mul(d); // (1+ε)² = 1+2ε (since ε²=0)
        let d_cubed = d_squared.mul(d); // (1+2ε)*(1+ε) = 1+3ε

        // 27a−d³ = (2+2ε) - (1+3ε) = 1-ε
        let term = twenty_seven_a.sub(d_cubed);

        // a(27a−d³) = (1+ε)*(1-ε) = 1-ε+ε-ε² = 1
        let condition = a.mul(term);

        // 1 is invertible in F5, so this condition is satisfied
        assert!(condition.is_invertible(), "a(27a−d³) must be invertible");

        // create the point P = [1, 2, 3+ε] from the paper
        let x = RingElement::from_field(field_1); // 1
        let y = RingElement::from_field(field_2); // 2
        let z = RingElement::new(field_3, field_1); // 3+ε

        let p = Projective::new(x, y, z);

        // verify P is on the curve
        assert!(p.is_on_curve(a, d), "P should be on the curve");

        // according to the paper, 4P = [1, 4, 3+2ε]
        let four_p = p.scalar_mul(4, a);

        let expected_4p_x = RingElement::from_field(field_1); // 1
        let expected_4p_y = RingElement::from_field(field_4); // 4
        let expected_4p_z = RingElement::new(field_3, field_2); // 3+2ε

        let expected_4p = Projective::new(expected_4p_x, expected_4p_y, expected_4p_z);

        assert!(
            four_p.is_equal(&expected_4p),
            "4P should equal [1, 4, 3+2ε]"
        );

        // according to the paper, 5P = [1, 3+2ε, 4+3ε]
        let five_p = p.scalar_mul(5, a);

        let expected_5p_x = RingElement::from_field(field_1); // 1
        let expected_5p_y = RingElement::new(field_3, field_2); // 3+2ε
        let expected_5p_z = RingElement::new(field_4, field_3); // 4+3ε

        let expected_5p = Projective::new(expected_5p_x, expected_5p_y, expected_5p_z);

        assert!(
            five_p.is_equal(&expected_5p),
            "5P should equal [1, 3+2ε, 4+3ε]"
        );

        // according to the paper, 35P = [1, 3, 2]
        let thirtyfive_p = p.scalar_mul(35, a);

        let expected_35p_x = RingElement::from_field(field_1); // 1
        let expected_35p_y = RingElement::from_field(field_3); // 3
        let expected_35p_z = RingElement::from_field(field_2); // 2

        let expected_35p =
            Projective::new(expected_35p_x, expected_35p_y, expected_35p_z);

        assert!(
            thirtyfive_p.is_equal(&expected_35p),
            "35P should equal [1, 3, 2]"
        );
    }

    #[test]
    fn kats_paper_3_2() {
        // using F11[ε] and the curve with a = 1+2ε, d = 2+ε as in Section 3.2
        type F11 = Fq<11>;

        // Create field elements
        let field_1 = F11::new(1);
        let field_2 = F11::new(2);
        let field_4 = F11::new(4);
        let field_6 = F11::new(6);
        let field_7 = F11::new(7);

        // parameters a = 1+2ε, d = 2+ε
        let a = RingElement::new(field_1, field_2); // 1+2ε
        let d = RingElement::new(field_2, field_1); // 2+ε

        assert!(
            Projective::verify_curve_constraints(a, d),
            "a(27a−d³) must be invertible"
        );

        // point P = [1, 7+6ε, 4+6ε] from the paper
        let x = RingElement::from_field(field_1); // 1
        let y = RingElement::new(field_7, field_6); // 7+6ε
        let z = RingElement::new(field_4, field_6); // 4+6ε

        let p = Projective::new(x, y, z);

        // verify P is on the curve
        assert!(p.is_on_curve(a, d), "P should be on the curve");
    }
}
