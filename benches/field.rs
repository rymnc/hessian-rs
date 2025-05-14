use core::ops::{
    Add,
    Mul,
    Sub,
};
use hessian_rs::Fq;
use rand::{
    Rng,
    thread_rng,
};

fn main() {
    divan::main();
}

#[crabtime::function]
fn bench_fq(moduli: Vec<u64>) {
    for modulus in moduli {
        crabtime::output! {
            fn generate_fq_{{modulus}}() -> Fq<{{modulus}}> {
                let mut rng = thread_rng();
                let value = rng.gen_range(1..i64::MAX as u64);
                Fq::<{{modulus}}>::new(value)
            }

            #[divan::bench]
            fn add_fq_{{modulus}}(bencher: divan::Bencher) {
                let mut rng = thread_rng();
                let augend = generate_fq_{{modulus}}();
                let addend = generate_fq_{{modulus}}();

                bencher.bench(|| {
                    augend.add(addend)
                });
            }

            #[divan::bench]
            fn sub_fq_{{modulus}}(bencher: divan::Bencher) {
                let mut rng = thread_rng();
                let minuend = generate_fq_{{modulus}}();
                let subtrahend = generate_fq_{{modulus}}();

                bencher.bench(|| {
                    minuend.sub(subtrahend)
                });
            }

            #[divan::bench]
            fn mul_fq_{{modulus}}(bencher: divan::Bencher) {
                let mut rng = thread_rng();
                let multiplicand = generate_fq_{{modulus}}();
                let multiplier = generate_fq_{{modulus}}();

                bencher.bench(|| {
                    multiplicand.mul(multiplier)
                });
            }

            #[divan::bench]
            fn inv_fq_{{modulus}}(bencher: divan::Bencher) {
                let mut rng = thread_rng();
                let element = generate_fq_{{modulus}}();

                bencher.bench(|| {
                    element.inv()
                });
            }

            #[divan::bench]
            fn pow_fq_{{modulus}}(bencher: divan::Bencher) {
                let mut rng = thread_rng();
                let element = generate_fq_{{modulus}}();
                let exponent = rng.gen_range(1..i64::MAX as u64);

                bencher.bench(|| {
                    element.pow(exponent)
                });
            }
        }
    }
}

// random selection of prime numbers
bench_fq!([41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97]);
