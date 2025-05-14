use core::ops::{
    Add,
    Mul,
    Sub,
};
use hessian_rs::{
    Fq,
    RingElement,
};
use rand::{
    Rng,
    thread_rng,
};

fn main() {
    divan::main();
}

#[crabtime::function]
fn bench_ring_element(moduli: Vec<u64>) {
    for modulus in moduli {
        crabtime::output! {
            fn generate_ring_element_{{modulus}}() -> RingElement<{{modulus}}> {
                let mut rng = thread_rng();
                let a_value = rng.gen_range(1..{{modulus}});
                let b_value = rng.gen_range(1..{{modulus}});

                RingElement::new(
                    Fq::<{{modulus}}>::new(a_value),
                    Fq::<{{modulus}}>::new(b_value)
                )
            }

            fn generate_field_element_{{modulus}}() -> Fq<{{modulus}}> {
                let mut rng = thread_rng();
                let value = rng.gen_range(1..{{modulus}});
                Fq::<{{modulus}}>::new(value)
            }

            #[divan::bench]
            fn add_ring_element_{{modulus}}(bencher: divan::Bencher) {
                let r1 = generate_ring_element_{{modulus}}();
                let r2 = generate_ring_element_{{modulus}}();

                bencher.bench(|| {
                    r1.add(r2)
                });
            }

            #[divan::bench]
            fn sub_ring_element_{{modulus}}(bencher: divan::Bencher) {
                let r1 = generate_ring_element_{{modulus}}();
                let r2 = generate_ring_element_{{modulus}}();

                bencher.bench(|| {
                    r1.sub(r2)
                });
            }

            #[divan::bench]
            fn mul_ring_element_{{modulus}}(bencher: divan::Bencher) {
                let r1 = generate_ring_element_{{modulus}}();
                let r2 = generate_ring_element_{{modulus}}();

                bencher.bench(|| {
                    r1.mul(r2)
                });
            }

            #[divan::bench]
            fn inv_ring_element_{{modulus}}(bencher: divan::Bencher) {
                let mut rng = thread_rng();
                let a_value = rng.gen_range(1..{{modulus}});
                let b_value = rng.gen_range(0..{{modulus}});

                let element = RingElement::new(
                    Fq::<{{modulus}}>::new(a_value),
                    Fq::<{{modulus}}>::new(b_value)
                );

                bencher.bench(|| {
                    element.inv()
                });
            }

            #[divan::bench]
            fn pow_ring_element_{{modulus}}(bencher: divan::Bencher) {
                let element = generate_ring_element_{{modulus}}();
                let mut rng = thread_rng();
                let exponent = rng.gen_range(1..100);

                bencher.bench(|| {
                    element.pow(exponent)
                });
            }

            #[divan::bench]
            fn from_field_ring_element_{{modulus}}(bencher: divan::Bencher) {
                let field_element = generate_field_element_{{modulus}}();

                bencher.bench(|| {
                    RingElement::from_field(field_element)
                });
            }

            #[divan::bench]
            fn is_invertible_ring_element_{{modulus}}(bencher: divan::Bencher) {
                let element = generate_ring_element_{{modulus}}();

                bencher.bench(|| {
                    element.is_invertible()
                });
            }


            #[divan::bench]
            fn equality_ring_element_{{modulus}}(bencher: divan::Bencher) {
                let r1 = generate_ring_element_{{modulus}}();
                let r2 = generate_ring_element_{{modulus}}();

                bencher.bench(|| {
                    r1 == r2
                });
            }
        }
    }
}

// random selection of prime numbers
bench_ring_element!([41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97]);
