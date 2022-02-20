//! Provides an implementation for jets with some utilities.
//!
//! A jet is n-dimensional dual number that can be used for automatic derivation.
//! See the [ceres library] for more information.
//!
//! [ceres library]: http://ceres-solver.org/automatic_derivatives.html#dual-numbers-jets
//!
//! # Feature flags
//!
//! These feature flags can be used to customize the crate.
//!
//! ## `big-rational` (disabled by default)
//!
//! If enabled, the type `num-rational::BigRational` will implement [`Number`].
//! This is very useful for writing tests with arbitrary precision. It will
//! require the dependency `num-rational`.

mod number;
pub use number::Number;

mod infinitesimal;
pub use infinitesimal::{Infinitesimal, NoInfinitesimal};

mod jet;
pub use jet::Jet;
