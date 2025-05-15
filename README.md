# hessian-rs

A Rust implementation of Twisted Hessian Curve cryptography over the local ring Fq[ε] where ε² = 0, with zero external runtime dependencies.

## Features

- Implementation of finite field arithmetic over Fq
- Local ring Fq[ε] implementation with ε² = 0
- Twisted Hessian curve operations in projective coordinates
- Diffie-Hellman key exchange protocol
- `no_std` compatible
- Zero dependencies for the core library
- Comprehensive test suite with known-answer tests from academic papers
- Benchmarking suite using crabtime and divan

## Implementation Details

This library implements twisted Hessian curves of the form:

$$
aX³ + Y³ + Z³ = dXYZ
$$

Where operations are performed over the local ring Fq[ε] with ε² = 0. Elements in this ring take the form a + bε where a, b ∈ Fq.

The implementation follows the mathematical foundations described in "Cryptography Over Twisted Hessian Curves of the Ring Fq[ε]" by Grini, Chillali, and Mouanis (2021).

## TODOs

The following areas are marked for improvement:

- Optimize finite field inversions using extended GCD
- Optimize power functions using Fermat's Little Theorem
- Optimize scalar multiplication using multi-scalar multiplication techniques
- Generate curve parameters for different finite fields (mentioned in `benches/curve.rs`)
- Additional testing and edge case handling

## Usage

Add this to your `Cargo.toml`:

```toml
[dependencies]
hessian-rs = "0.0.1"
```

### Example: Creating a Twisted Hessian curve

```rust
use hessian_rs::{Fq, RingElement, TwistedHessianCurve};

// Define a field with modulus 11
type F11 = Fq<11>;

// Create curve parameters a = 1+2ε, d = 2+ε
let a = RingElement::new(F11::new(1), F11::new(2)); // 1+2ε
let d = RingElement::new(F11::new(2), F11::new(1)); // 2+ε

// Create the curve
let curve = TwistedHessianCurve::new(a, d);
```

### Example: Diffie-Hellman Key Exchange

```rust
use hessian_rs::{Fq, RingElement, TwistedHessianCurve, Projective, dh::DiffieHellman};

// Define a field with modulus 5
type F5 = Fq<5>;

// Create curve parameters and a generator point
let a = RingElement::new(F5::new(1), F5::new(1)); // 1+ε
let d = RingElement::new(F5::new(1), F5::new(1)); // 1+ε
let curve = TwistedHessianCurve::new(a, d);

// Create the generator point P = [1, 2, 3+ε]
let generator = Projective::new(
    RingElement::from_field(F5::new(1)),
    RingElement::from_field(F5::new(2)),
    RingElement::new(F5::new(3), F5::new(1)),
);

// Create DH instance (order 45 as per the paper)
let dh = DiffieHellman::new(curve, generator, 45);

// Alice and Bob's private keys
let alice_private = 4;
let bob_private = 35;

// Generate key pairs
let (_, alice_public) = dh.generate_keypair(alice_private);
let (_, bob_public) = dh.generate_keypair(bob_private);

// Compute shared secrets
let alice_shared = dh.compute_shared_secret(alice_private, &bob_public);
let bob_shared = dh.compute_shared_secret(bob_private, &alice_public);

// alice_shared and bob_shared should be equal
assert!(alice_shared.is_equal(&bob_shared));
```

## Testing and Benchmarking

### Running Tests

```bash
cargo test
```

### Running Benchmarks

```bash
# Run all benchmarks
cargo bench

# Run specific benchmark
cargo bench --bench field
cargo bench --bench ring
cargo bench --bench projective
```

The benchmarks use [crabtime](https://github.com/khonsulabs/crabtime) for macro generation and [divan](https://github.com/nvzqz/divan) for measurement.

## References

- Grini, A., Chillali, A., & Mouanis, H. (2021). Cryptography over twisted Hessian curves of the ring Fq[ε]. Advances in Mathematics: Scientific Journal, 10(1), 235-243.
