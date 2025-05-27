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
    fn add__identity_element() {
        type F13 = Fq<13>;
        type R13 = RingElement<13>;

        let zero = R13::from_field(F13::new(0));
        let element = R13::new(F13::new(7), F13::new(4)); // 7 + 4ε

        // 0 + x = x
        let sum1 = zero.add(element);
        assert_eq!(sum1.constant().value(), 7);
        assert_eq!(sum1.epsilon_coeff().value(), 4);

        // x + 0 = x
        let sum2 = element.add(zero);
        assert_eq!(sum2.constant().value(), 7);
        assert_eq!(sum2.epsilon_coeff().value(), 4);
    }

    #[test]
    fn add__commutativity() {
        type F17 = Fq<17>;
        type R17 = RingElement<17>;

        let r1 = R17::new(F17::new(13), F17::new(8));
        let r2 = R17::new(F17::new(9), F17::new(15));

        let sum1 = r1.add(r2);
        let sum2 = r2.add(r1);

        assert_eq!(sum1.constant().value(), sum2.constant().value());
        assert_eq!(sum1.epsilon_coeff().value(), sum2.epsilon_coeff().value());
    }

    #[test]
    fn add__associativity() {
        type F19 = Fq<19>;
        type R19 = RingElement<19>;

        let r1 = R19::new(F19::new(5), F19::new(7));
        let r2 = R19::new(F19::new(11), F19::new(3));
        let r3 = R19::new(F19::new(17), F19::new(13));

        // (r1 + r2) + r3
        let sum1 = r1.add(r2).add(r3);
        // r1 + (r2 + r3)
        let sum2 = r1.add(r2.add(r3));

        assert_eq!(sum1.constant().value(), sum2.constant().value());
        assert_eq!(sum1.epsilon_coeff().value(), sum2.epsilon_coeff().value());
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
    fn sub__inverse_property() {
        type F23 = Fq<23>;
        type R23 = RingElement<23>;

        let r = R23::new(F23::new(15), F23::new(19));
        let zero = R23::from_field(F23::new(0));

        // r - r = 0
        let diff = r.sub(r);
        assert_eq!(diff.constant().value(), 0);
        assert_eq!(diff.epsilon_coeff().value(), 0);

        // 0 - r = -r
        let neg_r = zero.sub(r);
        // -15 ≡ 8 (mod 23), -19 ≡ 4 (mod 23)
        assert_eq!(neg_r.constant().value(), 8);
        assert_eq!(neg_r.epsilon_coeff().value(), 4);

        // r + (-r) = 0
        let sum = r.add(neg_r);
        assert_eq!(sum.constant().value(), 0);
        assert_eq!(sum.epsilon_coeff().value(), 0);
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
    fn mul__epsilon_squared_is_zero() {
        type F29 = Fq<29>;
        type R29 = RingElement<29>;

        let epsilon = R29::new(F29::new(0), F29::new(1)); // ε

        // ε * ε = 0
        let epsilon_squared = epsilon.mul(epsilon);
        assert_eq!(epsilon_squared.constant().value(), 0);
        assert_eq!(epsilon_squared.epsilon_coeff().value(), 0);

        // test with a more complex element
        let r = R29::new(F29::new(7), F29::new(13)); // 7 + 13ε
        let r_epsilon = r.mul(epsilon); // (7 + 13ε)ε = 7ε (since ε² = 0)
        assert_eq!(r_epsilon.constant().value(), 0);
        assert_eq!(r_epsilon.epsilon_coeff().value(), 7);
    }

    #[test]
    fn mul__zero_element() {
        type F37 = Fq<37>;
        type R37 = RingElement<37>;

        let zero = R37::from_field(F37::new(0));
        let r = R37::new(F37::new(25), F37::new(31));

        // 0 * r = 0
        let prod1 = zero.mul(r);
        assert_eq!(prod1.constant().value(), 0);
        assert_eq!(prod1.epsilon_coeff().value(), 0);

        // r * 0 = 0
        let prod2 = r.mul(zero);
        assert_eq!(prod2.constant().value(), 0);
        assert_eq!(prod2.epsilon_coeff().value(), 0);
    }

    #[test]
    fn mul__distributivity() {
        type F41 = Fq<41>;
        type R41 = RingElement<41>;

        let r1 = R41::new(F41::new(7), F41::new(13));
        let r2 = R41::new(F41::new(19), F41::new(23));
        let r3 = R41::new(F41::new(31), F41::new(37));

        // r1 * (r2 + r3) = r1 * r2 + r1 * r3
        let left = r1.mul(r2.add(r3));
        let right = r1.mul(r2).add(r1.mul(r3));

        assert_eq!(left.constant().value(), right.constant().value());
        assert_eq!(left.epsilon_coeff().value(), right.epsilon_coeff().value());
    }

    #[test]
    fn mul__associativity() {
        type F43 = Fq<43>;
        type R43 = RingElement<43>;

        let r1 = R43::new(F43::new(11), F43::new(17));
        let r2 = R43::new(F43::new(23), F43::new(29));
        let r3 = R43::new(F43::new(37), F43::new(41));

        // (r1 * r2) * r3
        let prod1 = r1.mul(r2).mul(r3);
        // r1 * (r2 * r3)
        let prod2 = r1.mul(r2.mul(r3));

        assert_eq!(prod1.constant().value(), prod2.constant().value());
        assert_eq!(prod1.epsilon_coeff().value(), prod2.epsilon_coeff().value());
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

    #[test]
    #[should_panic(expected = "Element not invertible")]
    fn inv__non_invertible_element() {
        type F53 = Fq<53>;
        type R53 = RingElement<53>;

        // Element with a = 0 is not invertible
        let r = R53::new(F53::new(0), F53::new(25));
        r.inv(); // Should panic
    }

    #[test]
    fn complex_arithmetic_chains() {
        type F71 = Fq<71>;
        type R71 = RingElement<71>;

        let r1 = R71::new(F71::new(15), F71::new(23));
        let r2 = R71::new(F71::new(37), F71::new(41));
        let r3 = R71::new(F71::new(53), F71::new(59));

        // (r1 + r2) * r3 - r1 * (r2 - r3)
        let expr1 = r1.add(r2).mul(r3);
        let expr2 = r1.mul(r2.sub(r3));
        let result = expr1.sub(expr2);

        // r1*r3 + r2*r3 - r1*r2 + r1*r3 = 2*r1*r3 + r2*r3 - r1*r2
        let verify = r1.mul(r3).add(r1.mul(r3)).add(r2.mul(r3)).sub(r1.mul(r2));

        assert_eq!(result.constant().value(), verify.constant().value());
        assert_eq!(
            result.epsilon_coeff().value(),
            verify.epsilon_coeff().value()
        );
    }
}
