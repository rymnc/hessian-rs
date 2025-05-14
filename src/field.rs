//! Finite field implementation

use core::ops::{
    Add,
    Mul,
    Sub,
};

/// Finite field Fq implementation where q is a prime power
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Fq<const Q: u64> {
    value: u64,
}

impl<const Q: u64> Fq<Q> {
    /// Create a new element in the finite field Fq
    pub const fn new(value: u64) -> Self {
        assert!(Q > 1, "Field modulus must be greater than 1");
        assert!(
            Q < i64::MAX as u64,
            "Field modulus must be less than i64::MAX"
        );
        assert!(value < i64::MAX as u64, "value must be less than i64::MAX");

        let value = value.rem_euclid(Q);
        Fq { value }
    }

    /// Get the value of the field element
    pub fn value(&self) -> u64 {
        self.value
    }

    /// Get the modulus of the field
    pub fn modulus() -> u64 {
        Q
    }

    /// Multiplicative inverse of a field element
    pub fn inv(&self) -> Self {
        // TODO: optimize using extended gcd
        assert_ne!(self.value, 0, "Cannot invert zero");

        let mut s = 0i64;
        let mut old_s = 1i64;
        let mut t = 1i64;
        let mut old_t = 0i64;
        let mut r = Q as i64;
        let mut old_r = self.value as i64;

        while r != 0 {
            let quotient = old_r.checked_div(r).expect("division failed");

            let temp = old_r;
            old_r = r;
            r = temp
                .checked_sub(quotient.checked_mul(r).expect("multiplication failed"))
                .expect("subtraction failed");

            let temp = old_s;
            old_s = s;
            s = temp
                .checked_sub(quotient.checked_mul(s).expect("multiplication failed"))
                .expect("subtraction failed");

            let temp = old_t;
            old_t = t;
            t = temp
                .checked_sub(quotient.checked_mul(t).expect("multiplication failed"))
                .expect("subtraction failed");
        }

        // if old_r > 1, then gcd(a, m) != 1 and inverse doesn't exist
        assert_eq!(old_r, 1, "Inverse doesn't exist (gcd != 1)");

        let result = if old_s < 0 {
            old_s
                .checked_add(i64::try_from(Q).unwrap()) // qed
                .expect("addition failed")
        } else {
            old_s
        };
        Fq::new(result as u64)
    }

    /// Pow
    pub fn pow(&self, exponent: u64) -> Self {
        // TODO: optimize using fermat's little theorem
        let mut result = Fq::new(1);
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

    /// This is needed for twisted Hessian curve conditions
    pub fn is_minus_three_square() -> bool {
        // -3 mod p = p-3 mod p
        let minus_three = Fq::<Q>::new(Q.checked_sub(3).expect("subtraction failed"));

        // a^((p-1)/2) ≡ 1 mod p, if a is a quadratic residue
        if Q % 2 == 0 {
            return false;
        }

        let exponent = (Q.checked_sub(1).expect("subtraction failed"))
            .checked_div(2)
            .expect("division failed");
        minus_three.pow(exponent).value == 1
    }
}

impl<const Q: u64> Add for Fq<Q> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let sum =
            (self.value.checked_add(rhs.value).expect("addition failed")).rem_euclid(Q);
        Fq::new(sum)
    }
}

#[allow(clippy::arithmetic_side_effects)]
impl<const Q: u64> Add for &Fq<Q> {
    type Output = Fq<Q>;

    fn add(self, rhs: Self) -> Self::Output {
        *self + *rhs
    }
}

#[allow(clippy::arithmetic_side_effects)]
impl<const Q: u64> Sub for Fq<Q> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        let mut diff = self.value as i64 - rhs.value as i64;
        if diff < 0 {
            diff += Q as i64;
        }
        Fq::new(diff as u64)
    }
}

#[allow(clippy::arithmetic_side_effects)]
impl<const Q: u64> Sub for &Fq<Q> {
    type Output = Fq<Q>;

    fn sub(self, rhs: Self) -> Self::Output {
        *self - *rhs
    }
}

impl<const Q: u64> Mul for Fq<Q> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let product = (self
            .value
            .checked_mul(rhs.value)
            .expect("multiplication failed"))
        .rem_euclid(Q);
        Fq::new(product)
    }
}

#[allow(clippy::arithmetic_side_effects)]
impl<const Q: u64> Mul for &Fq<Q> {
    type Output = Fq<Q>;

    fn mul(self, rhs: Self) -> Self::Output {
        *self * *rhs
    }
}

#[allow(non_snake_case)]
#[cfg(test)]
mod tests {
    use proptest::proptest;

    use super::*;

    #[test]
    fn add__computes_correctly_without_overflow() {
        type F11 = Fq<11>;

        let a = F11::new(5);
        let b = F11::new(3);

        assert_eq!(a.add(b).value(), 8); // 5 + 3 = 8
    }

    #[test]
    fn add__computes_correctly_with_overflow() {
        type F11 = Fq<11>;

        let a = F11::new(7);
        let b = F11::new(8);
        assert_eq!(a.add(b).value(), 4); // 7 + 8 = 15 ≡ 4 (mod 11)
    }

    #[test]
    fn sub__computes_correctly_without_underflow() {
        type F11 = Fq<11>;

        let a = F11::new(5);
        let b = F11::new(3);

        assert_eq!(a.sub(b).value(), 2); // 5 - 3 = 2
    }

    #[test]
    fn sub__computes_correctly_with_underflow() {
        type F11 = Fq<11>;

        let a = F11::new(2);
        let b = F11::new(5);
        assert_eq!(a.sub(b).value(), 8); // 2 - 5 = -3 ≡ 8 (mod 11)
    }

    #[test]
    fn mul__computes_correctly_without_overflow() {
        type F11 = Fq<11>;

        let a = F11::new(5);
        let b = F11::new(3);

        assert_eq!(a.mul(b).value(), 4); // 5 * 3 = 15 ≡ 4 (mod 11)
    }

    #[test]
    fn mul__computes_correctly_with_overflow() {
        type F11 = Fq<11>;

        let a = F11::new(6);
        let b = F11::new(9);
        assert_eq!(a.mul(b).value(), 10); // 6 * 9 = 54 ≡ 10 (mod 11)
    }

    #[test]
    fn inv__computes_correctly() {
        type F11 = Fq<11>;

        let a = F11::new(5);
        let a_inv = a.inv();

        assert_eq!(a.mul(a_inv).value(), 1);
    }

    #[test]
    fn inv__F11_kats() {
        type F11 = Fq<11>;

        assert_eq!(F11::new(1).inv().value(), 1);
        assert_eq!(F11::new(2).inv().value(), 6);
        assert_eq!(F11::new(3).inv().value(), 4);
        assert_eq!(F11::new(4).inv().value(), 3);
        assert_eq!(F11::new(5).inv().value(), 9);
        assert_eq!(F11::new(6).inv().value(), 2);
    }

    #[test]
    fn inv__proptest() {
        // has to be prime
        const MAX_FIELD: u64 = 7919;
        type ProptestField = Fq<MAX_FIELD>;

        proptest!(|(a in 1..MAX_FIELD)| {
            let a = ProptestField::new(a);
            let a_inv = a.inv();
            assert_eq!(a.mul(a_inv).value(), 1);
        });
    }

    #[test]
    #[should_panic(expected = "Cannot invert zero")]
    fn inv__should_panic_when_inverting_zero() {
        let zero = Fq::<11>::new(0);
        zero.inv();
    }

    #[test]
    #[should_panic(expected = "Inverse doesn't exist (gcd != 1)")]
    fn inv__should_panic_when_mod_is_not_prime() {
        let two = Fq::<4>::new(2);
        assert_eq!(two.inv().value(), 2);
    }

    #[test]
    fn pow__computes_correctly_without_overflow() {
        type F11 = Fq<11>;

        let a = F11::new(2);

        // 2^0 = 1
        assert_eq!(a.pow(0).value(), 1);

        // 2^1 = 2
        assert_eq!(a.pow(1).value(), 2);

        // 2^3 = 8
        assert_eq!(a.pow(3).value(), 8);
    }

    #[test]
    fn pow__computes_correctly_with_overflow() {
        type F11 = Fq<11>;

        let a = F11::new(2);

        // 2^10 = 1024 ≡ 1 (mod 11)
        assert_eq!(a.pow(10).value(), 1);
    }

    #[test]
    fn minus_three_square_in_field5() {
        // In F5, -3 ≡ 2 (mod 5)
        // 2^((5-1)/2) = 2^2 = 4 ≡ -1 (mod 5)
        // So 2 is not a quadratic residue in F5
        assert!(!Fq::<5>::is_minus_three_square());
    }

    #[test]
    fn minus_three_square_in_field7() {
        // In F7, -3 ≡ 4 (mod 7)
        // 4 = 2^2, so it's a perfect square in F7
        assert!(Fq::<7>::is_minus_three_square());
    }

    #[test]
    fn minus_three_square_in_various_fields() {
        assert!(!Fq::<11>::is_minus_three_square()); // -3 ≡ 8 (mod 11)
        assert!(Fq::<13>::is_minus_three_square()); // -3 ≡ 10 (mod 13)
        assert!(!Fq::<17>::is_minus_three_square()); // -3 ≡ 14 (mod 17)
        assert!(Fq::<19>::is_minus_three_square()); // -3 ≡ 16 (mod 19)
    }
}
