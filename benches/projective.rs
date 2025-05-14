use hessian_rs::{
    Fq,
    Projective,
    RingElement,
    TwistedHessianCurve,
};
use rand::{
    Rng,
    thread_rng,
};

fn main() {
    divan::main();
}

#[crabtime::function]
fn bench_projective(moduli: Vec<u64>) {
    for modulus in moduli {
        crabtime::output! {
            fn create_curve_{{modulus}}() -> TwistedHessianCurve<{{modulus}}> {
                let field_1 = Fq::<{{modulus}}>::new(1);
                let field_2 = Fq::<{{modulus}}>::new(2);

                let a = RingElement::new(field_1, field_2);  // a = 1+2ε
                let d = RingElement::new(field_2, field_1);  // d = 2+ε

                TwistedHessianCurve::new(a, d)
            }

            fn generate_point_{{modulus}}() -> Projective<{{modulus}}> {
                let curve = create_curve_{{modulus}}();

                let mut rng = thread_rng();

                for _ in 0..100 {
                    let x_val = rng.gen_range(1..{{modulus}});
                    let y_val = rng.gen_range(1..{{modulus}});
                    let z_val = rng.gen_range(1..{{modulus}});

                    let x = RingElement::from_field(Fq::<{{modulus}}>::new(x_val));
                    let y = RingElement::from_field(Fq::<{{modulus}}>::new(y_val));
                    let z = RingElement::from_field(Fq::<{{modulus}}>::new(z_val));

                    let point = Projective::new(x, y, z);

                    if curve.contains(&point) {
                        return point;
                    }
                }

                let field_1 = Fq::<{{modulus}}>::new(1);
                Projective::new(
                    RingElement::from_field(field_1),
                    RingElement::from_field(field_1),
                    RingElement::from_field(field_1)
                )
            }

            fn generate_curve_parameter_{{modulus}}() -> RingElement<{{modulus}}> {
                let curve = create_curve_{{modulus}}();
                curve.a()
            }

            #[divan::bench]
            fn add_projective_{{modulus}}(bencher: divan::Bencher) {
                let p1 = generate_point_{{modulus}}();
                let p2 = generate_point_{{modulus}}();
                let a = generate_curve_parameter_{{modulus}}();

                bencher.bench(|| {
                    p1.add(&p2, a)
                });
            }

            #[divan::bench]
            fn double_projective_{{modulus}}(bencher: divan::Bencher) {
                let p = generate_point_{{modulus}}();
                let a = generate_curve_parameter_{{modulus}}();

                bencher.bench(|| {
                    p.double(a)
                });
            }

            #[divan::bench]
            fn scalar_mul_projective_{{modulus}}(bencher: divan::Bencher) {
                let p = generate_point_{{modulus}}();
                let a = generate_curve_parameter_{{modulus}}();
                let scalar = 2; // TODO: reduce flake / randomness

                bencher.bench(|| {
                    p.scalar_mul(scalar, a)
                });
            }

            #[divan::bench]
            fn negate_projective_{{modulus}}(bencher: divan::Bencher) {
                let p = generate_point_{{modulus}}();

                bencher.bench(|| {
                    p.negate()
                });
            }

            #[divan::bench]
            fn is_on_curve_projective_{{modulus}}(bencher: divan::Bencher) {
                let p = generate_point_{{modulus}}();
                let curve = create_curve_{{modulus}}();

                bencher.bench(|| {
                    curve.contains(&p)
                });
            }

            #[divan::bench]
            fn is_equal_projective_{{modulus}}(bencher: divan::Bencher) {
                let p1 = generate_point_{{modulus}}();
                let p2 = generate_point_{{modulus}}();

                bencher.bench(|| {
                    p1.is_equal(&p2)
                });
            }
        }
    }
}

// random selection of prime numbers
bench_projective!([41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97]);
