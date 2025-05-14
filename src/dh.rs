//! ECDH
use crate::{
    curve::TwistedHessianCurve,
    projective::Projective,
};

/// ECDH for a twisted hessian curve w/ ring
pub struct DiffieHellman<const Q: u64> {
    curve: TwistedHessianCurve<Q>,
    generator: Projective<Q>,
    order: u64,
}

impl<const Q: u64> DiffieHellman<Q> {
    /// New ECDH with provided generator point and curve
    pub fn new(
        curve: TwistedHessianCurve<Q>,
        generator: Projective<Q>,
        order: u64,
    ) -> Self {
        assert!(curve.contains(&generator), "Generator must be on the curve");

        // Verify the order is correct
        let identity = curve.identity();
        let check = curve.scalar_mul(&generator, order);
        assert!(
            check.is_equal(&identity),
            "Generator's order must match the provided order"
        );

        DiffieHellman {
            curve,
            generator,
            order,
        }
    }

    /// Generate a new key pair (private key, public key)
    pub fn generate_keypair(&self, private_key: u64) -> (u64, Projective<Q>) {
        // Ensure private key is within the valid range
        let private_key = private_key % self.order;
        if private_key == 0 {
            panic!("Private key cannot be zero");
        }

        let public_key = self.curve.scalar_mul(&self.generator, private_key);
        (private_key, public_key)
    }

    /// Compute the shared secret from a private key and another party's public key
    pub fn compute_shared_secret(
        &self,
        private_key: u64,
        public_key: &Projective<Q>,
    ) -> Projective<Q> {
        assert!(
            self.curve.contains(public_key),
            "Public key must be on the curve"
        );
        self.curve.scalar_mul(public_key, private_key)
    }
}

/// Simulates a Diffie-Hellman key exchange between two parties
pub fn simulate_key_exchange<const Q: u64>(
    dh: &DiffieHellman<Q>,
    alice_private: u64,
    bob_private: u64,
) -> (Projective<Q>, Projective<Q>) {
    let (_, alice_public) = dh.generate_keypair(alice_private);
    let (_, bob_public) = dh.generate_keypair(bob_private);

    let alice_shared = dh.compute_shared_secret(alice_private, &bob_public);
    let bob_shared = dh.compute_shared_secret(bob_private, &alice_public);

    (alice_shared, bob_shared)
}

#[allow(non_snake_case)]
#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        field::Fq,
        ring::RingElement,
    };

    #[test]
    fn diffie_hellman__should_succeed__kats() {
        // Section 3.1.2 example: Using F5[ε] with a = 1+ε, d = 1+ε
        type F5 = Fq<5>;

        let field_1 = F5::new(1);
        let field_2 = F5::new(2);
        let field_3 = F5::new(3);

        // create curve parameters a = 1+ε, d = 1+ε
        let a = RingElement::new(field_1, field_1); // 1+ε
        let d = RingElement::new(field_1, field_1); // 1+ε

        let curve = TwistedHessianCurve::new(a, d);

        // create the generator point P = [1, 2, 3+ε] from the paper
        let x = RingElement::from_field(field_1); // 1
        let y = RingElement::from_field(field_2); // 2
        let z = RingElement::new(field_3, field_1); // 3+ε

        let generator = Projective::new(x, y, z);

        assert!(
            curve.contains(&generator),
            "Generator point must be on the curve"
        );

        // the paper says this point has order 45
        let order = 45;

        let dh = DiffieHellman::new(curve, generator, order);

        let alice_private = 4;
        let bob_private = 35;

        // simulate key exchange
        let (alice_shared, bob_shared) =
            simulate_key_exchange(&dh, alice_private, bob_private);

        assert!(
            alice_shared.is_equal(&bob_shared),
            "Shared secrets should match"
        );

        // the expected shared secret is 5P = [1, 3+2ε, 4+3ε]
        let expected_secret = Projective::new(
            RingElement::from_field(field_1),      // 1
            RingElement::new(field_3, F5::new(2)), // 3+2ε
            RingElement::new(F5::new(4), field_3), // 4+3ε
        );

        assert!(
            alice_shared.is_equal(&expected_secret),
            "Shared secret should equal expected value from paper"
        );
    }
}
