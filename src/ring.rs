//! Ring implementation for Fq[ε] where ε² = 0

use crate::field::Fq;
use core::ops::Mul;

/// Element in the local ring Fq[ε] where ε² = 0
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RingElement<const Q: u64> {
    a: Fq<Q>,
    b: Fq<Q>,
}

impl<const Q: u64> RingElement<Q> {
    /// Create a new element a + bε in the local ring Fq[ε]
    pub const fn new(a: Fq<Q>, b: Fq<Q>) -> Self {
        RingElement { a, b }
    }

    /// Create an element a in the local ring (without ε component)
    pub const fn from_field(a: Fq<Q>) -> Self {
        RingElement::new(a, Fq::new(0))
    }

    /// Get the constant part (a) of a + bε
    pub fn constant(&self) -> Fq<Q> {
        self.a
    }

    /// Get the coefficient (b) of ε in a + bε
    pub fn epsilon_coeff(&self) -> Fq<Q> {
        self.b
    }

    /// Get the modulus of the underlying field
    pub fn modulus() -> u64 {
        Q
    }

    /// Check if this ring element is invertible
    pub fn is_invertible(&self) -> bool {
        // a + bε is invertible if a is non-zero in Fq
        self.a.value() != 0
    }

    /// Multiplicative inverse of a ring element
    pub fn inv(&self) -> Self {
        assert!(self.is_invertible(), "Element not invertible");

        // For a + bε, the inverse is a⁻¹ - ba⁻²ε
        let a_inv = self.a.inv();
        let a_inv_squared = a_inv.mul(a_inv);
        let b_a_inv_squared = self.b.mul(a_inv_squared);

        RingElement::new(
            a_inv,
            Fq::new(
                Q.checked_sub(b_a_inv_squared.value())
                    .expect("subtraction failed"),
            ),
        )
    }

    /// Raise a ring element to a power
    pub fn pow(&self, exponent: u64) -> Self {
        // TODO: optimize using fermat's little theorem
        let mut result = RingElement::new(Fq::new(1), Fq::new(0));
        let mut base = *self;
        let mut exp = exponent;

        while exp > 0 {
            if exp & 1 == 1 {
                result = result.mul(base);
            }
            base = base.mul(base);
            exp >>= 1;
        }

        result
    }
}

#[allow(clippy::arithmetic_side_effects)]
impl<const Q: u64> core::ops::Add for RingElement<Q> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        RingElement::new(self.a + other.a, self.b + other.b)
    }
}

#[allow(clippy::arithmetic_side_effects)]
impl<const Q: u64> core::ops::Sub for RingElement<Q> {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        RingElement::new(self.a - other.a, self.b - other.b)
    }
}

#[allow(clippy::arithmetic_side_effects)]
impl<const Q: u64> core::ops::Mul for RingElement<Q> {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        // (a + bε) * (c + dε) = ac + (ad + bc)ε + bdε² = ac + (ad + bc)ε
        // since ε² = 0
        let ac = self.a * other.a;
        let ad = self.a * other.b;
        let bc = self.b * other.a;

        RingElement::new(ac, ad + bc)
    }
}

#[allow(non_snake_case)]
#[cfg(test)]
mod tests {
    use super::*;
    use core::ops::{
        Add,
        Mul,
        Sub,
    };

    #[test]
    fn add__computes_correctly() {
        type F11 = Fq<11>;
        type R11 = RingElement<11>;

        let a_field = F11::new(5);
        let b_field = F11::new(3);
        let c_field = F11::new(2);
        let d_field = F11::new(7);

        let r1 = R11::new(a_field, b_field); // 5 + 3ε
        let r2 = R11::new(c_field, d_field); // 2 + 7ε

        // (5 + 3ε) + (2 + 7ε) = 7 + 10ε
        let r_sum = r1.add(r2);
        assert_eq!(r_sum.constant().value(), 7);
        assert_eq!(r_sum.epsilon_coeff().value(), 10);
    }

    #[test]
    fn sub__computes_correctly() {
        type F11 = Fq<11>;
        type R11 = RingElement<11>;

        let a_field = F11::new(5);
        let b_field = F11::new(3);
        let c_field = F11::new(2);
        let d_field = F11::new(7);

        let r1 = R11::new(a_field, b_field); // 5 + 3ε
        let r2 = R11::new(c_field, d_field); // 2 + 7ε

        // (5 + 3ε) - (2 + 7ε) = 3 - 4ε
        let r_diff = r1.sub(r2);
        assert_eq!(r_diff.constant().value(), 3);
        assert_eq!(r_diff.epsilon_coeff().value(), 7); // -4 ≡ 7 mod 11
    }

    #[test]
    fn mul__computes_correctly() {
        type F11 = Fq<11>;
        type R11 = RingElement<11>;

        let a_field = F11::new(5);
        let b_field = F11::new(3);
        let c_field = F11::new(2);
        let d_field = F11::new(7);

        let r1 = R11::new(a_field, b_field); // 5 + 3ε
        let r2 = R11::new(c_field, d_field); // 2 + 7ε

        // (5 + 3ε) * (2 + 7ε) = 10 + (5*7 + 3*2)ε
        let r_prod = r1.mul(r2);
        assert_eq!(r_prod.constant().value(), 10);
        assert_eq!(r_prod.epsilon_coeff().value(), 8); // 41 ≡ 8 mod 11
    }

    #[test]
    fn inv__computes_correctly() {
        type F11 = Fq<11>;
        type R11 = RingElement<11>;

        let a_field = F11::new(5);
        let b_field = F11::new(3);

        let r1 = R11::new(a_field, b_field); // 5 + 3ε

        // inverse of r1 = 5 + 3ε
        // for a + bε, the inverse is a⁻¹ - ba⁻²ε
        // 5⁻¹ ≡ 9 (mod 11)
        // -3*9² ≡ -3*81 ≡ -3*4 ≡ -12 ≡ 10 (mod 11)
        // inverse =  9 + 10ε
        let r1_inv = r1.inv();
        assert_eq!(r1_inv.constant().value(), 9);
        assert_eq!(r1_inv.epsilon_coeff().value(), 10);

        // verify r1 * r1_inv = 1 + 0ε
        let r_one = r1.mul(r1_inv);
        assert_eq!(r_one.constant().value(), 1);
        assert_eq!(r_one.epsilon_coeff().value(), 0);
    }
}
