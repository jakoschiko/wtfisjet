//! Provides an implementation for jets with some utilities.
//!
//! A jet is n-dimensional dual number that can be used for automatic derivation.
//! See the [ceres library] for more information.
//!
//! [ceres library]: http://ceres-solver.org/automatic_derivatives.html#dual-numbers-jets

mod number;
pub use number::Number;

mod infinitesimal;
pub use infinitesimal::Infinitesimal;

mod jet;
pub use jet::Jet;
