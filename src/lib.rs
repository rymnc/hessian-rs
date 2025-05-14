//! Twisted Hessian curve crypto over a ring

#![cfg_attr(not(test), no_std)]
#![deny(clippy::arithmetic_side_effects)]
#![deny(clippy::cast_possible_truncation)]
#![deny(unused_crate_dependencies)]
#![deny(missing_docs)]
#![deny(warnings)]

pub mod curve;
pub mod dh;
pub mod field;
pub mod projective;
pub mod ring;

// convenient re-exports
pub use curve::TwistedHessianCurve;
pub use field::Fq;
pub use projective::Projective;
pub use ring::RingElement;

#[cfg(test)]
use crabtime as _;
#[cfg(test)]
use divan as _;
#[cfg(test)]
use rand as _;
