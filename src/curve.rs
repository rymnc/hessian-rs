//! Twisted hessian curve over the ring Fq[ε]

use crate::{
    field::Fq,
    projective::Projective,
    ring::RingElement,
};
use core::ops::{
    Mul,
    Sub,
};

/// Represents a twisted Hessian curve aX³ + Y³ + Z³ = dXYZ over the ring Fq[ε]
#[derive(Debug, Clone, Copy)]
pub struct TwistedHessianCurve<const Q: u64> {
    a: RingElement<Q>,
    d: RingElement<Q>,
}

impl<const Q: u64> TwistedHessianCurve<Q> {
    /// Create a new twisted Hessian curve with parameters a and d
    pub fn new(a: RingElement<Q>, d: RingElement<Q>) -> Self {
        // check if a*(27a-d³) is invertible in R2
        let d_cubed = d.mul(d).mul(d);
        let twenty_seven = RingElement::from_field(Fq::new(27u64.rem_euclid(Q)));
        let term = twenty_seven.mul(a).sub(d_cubed);
        let check = a.mul(term);

        assert!(
            check.is_invertible(),
            "a*(27a-d³) must be invertible for a valid curve"
        );

        TwistedHessianCurve { a, d }
    }

    /// Get the a parameter of the curve
    pub fn a(&self) -> RingElement<Q> {
        self.a
    }

    /// Get the d parameter of the curve
    pub fn d(&self) -> RingElement<Q> {
        self.d
    }

    /// Get the modulus of the underlying field
    pub fn modulus(&self) -> u64 {
        Q
    }

    /// Get the identity element of the curve group
    pub fn identity(&self) -> Projective<Q> {
        Projective::identity()
    }

    /// Check if a point lies on this curve
    pub fn contains(&self, point: &Projective<Q>) -> bool {
        point.is_on_curve(self.a, self.d)
    }

    /// Add two points on this curve
    pub fn add(&self, p: &Projective<Q>, q: &Projective<Q>) -> Projective<Q> {
        assert!(self.contains(p), "Projective P must be on the curve");
        assert!(self.contains(q), "Projective Q must be on the curve");

        p.add(q, self.a)
    }

    /// Multiply a point by a scalar
    pub fn scalar_mul(&self, p: &Projective<Q>, scalar: u64) -> Projective<Q> {
        assert!(self.contains(p), "Projective must be on the curve");

        p.scalar_mul(scalar, self.a)
    }

    /// Calculate the order of a point (the smallest positive k such that k*P = O)
    pub fn point_order(&self, point: &Projective<Q>) -> u64 {
        // TODO: optimize this, rlc
        assert!(self.contains(point), "Projective must be on the curve");

        let identity = self.identity();

        // handle the case where the point is already the identity
        if point.is_equal(&identity) {
            return 1;
        }

        // for a curve over F_q[ε], the group order is <= q^2
        let max_possible_order = Q.pow(2);

        for order in 2..=max_possible_order {
            let multiple = self.scalar_mul(point, order);
            if multiple.is_equal(&identity) {
                let previous = self.scalar_mul(point, order - 1);
                if !previous.is_equal(&identity) {
                    return order;
                }
            }
        }

        // this should never happen for a valid point on the curve
        panic!(
            "Could not determine the order of the point within range 1..{}",
            max_possible_order
        );
    }
}

// TODO: more test cases
#[allow(non_snake_case)]
#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        field::Fq,
        ring::RingElement,
    };

    #[test]
    fn new__should_succeed__when__parameters_are_valid() {
        type F5 = Fq<5>;
        let field_1 = F5::new(1);
        let field_2 = F5::new(2);

        let a = RingElement::from_field(field_1);
        let d = RingElement::from_field(field_2);

        let curve = TwistedHessianCurve::new(a, d);
        assert_eq!(curve.a().constant().value(), 1);
        assert_eq!(curve.d().constant().value(), 2);
        assert_eq!(curve.modulus(), 5);
    }

    #[test]
    #[should_panic(expected = "a*(27a-d³) must be invertible for a valid curve")]
    fn new__should_fail__when__parameters_are_invalid() {
        // a=1, d=3 -> a(27a-d³) = 1(2-2) = 0 not invertible
        type F5 = Fq<5>;
        let a = RingElement::from_field(F5::new(1));
        let d_invalid = RingElement::from_field(F5::new(3));

        TwistedHessianCurve::new(a, d_invalid);
    }
}
