#![feature(generic_associated_types)]

//! Provides an implementation for jets with some utilities.
//!
//! A jet is n-dimensional dual number that can be used for automatic derivation.
//! See the [ceres library] for more information.
//!
//! [ceres library]: http://ceres-solver.org/automatic_derivatives.html#dual-numbers-jets
//!
//! # Nightly rust
//!
//! This crate uses the nightly feature `generic_associated_types`. Hence it's necessary to use
//! the nightly rust compiler. Hopefully this won't be necessary anymore in the near future,
//! see [The push for GATs stabilization].
//!
//! [The push for GATs stabilization]: https://blog.rust-lang.org/2021/08/03/GATs-stabilization-push.html
//!
//! # Feature flags
//!
//! These feature flags can be used to customize the crate.
//!
//! ## `big-rational-number` (disabled by default)
//!
//! If enabled, the type `num-rational::BigRational` will implement [`Number`].
//! This is very useful for writing tests with arbitrary precision. It will
//! require the dependency `num-rational`.

mod number;
pub use number::Number;

mod infinitesimal;
pub use infinitesimal::{DenseInfinitesimal, Infinitesimal, NoInfinitesimal};

mod jet;
pub use jet::Jet;
