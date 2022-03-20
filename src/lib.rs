#![feature(generic_associated_types)]
#![feature(total_cmp)]

//! Provides an implementation for jets with some utilities.
//!
//! A jet is n-dimensional dual number that can be used for automatic derivation.
//! See the [ceres library] for more information.
//!
//! [ceres library]: http://ceres-solver.org/automatic_derivatives.html#dual-numbers-jets
//!
//! # Nightly rust
//!
//! This crate uses some nightly features. Hence it's necessary to use the nightly rust compiler.
//! Hopefully this won't be necessary anymore in the future.
//!
//! Nightly features used by this crate:
//! - `generic_associated_types` is used for the inner [`Iterator`] types of [`Infinitesimal`].
//! We can expect the stabilization in the near future, see [The push for GATs stabilization].
//! - `total_cmp` allows us to use function like [`slice::binary_search_by`] for [`f32`] and
//! [`f64`]. Although it's handy, this could also be solved with workarounds. See the
//! [tracking issue for `total_cmp`].
//!
//! [The push for GATs stabilization]: https://blog.rust-lang.org/2021/08/03/GATs-stabilization-push.html
//! [tracking issue for `total_cmp`]: https://github.com/rust-lang/rust/issues/72599
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
//!
//! ## `sparse-infinitesimal` (disabled by default)
//!
//! If enabled, the type `wtfisjet::SparseInfinitesimal` is available. It will
//! require the dependency `intmap`.
//!
//! ## `const-infinitesimal` (disabled by default)
//!
//! If enabled, the type `wtfisjet::ConstInfinitesimal` is available. It will
//! require the dependency `array-init`.
//!
//! ## `dice` (disabled by default)
//!
//! If enabled, the module `wtfisjet::dice` is available. It provides generators for
//! random test data. It will require the dependency `dicetest`.
//!
//! ## `asserts` (disabled by default)
//!
//! If enabled, the module `wtfisjet::asserts` is available. It provides assertions for tests.
//! It will require the dependency `dicetest` and `diceprop`.

mod number;
pub use number::Number;

mod dim;
pub use dim::Dim;

mod infinitesimal;
#[cfg(any(test, feature = "const-infinitesimal"))]
pub use infinitesimal::ConstInfinitesimal;
#[cfg(any(test, feature = "sparse-infinitesimal"))]
pub use infinitesimal::SparseInfinitesimal;
pub use infinitesimal::{DenseInfinitesimal, Infinitesimal, NoInfinitesimal};

mod jet;
pub use jet::Jet;

#[cfg(any(test, feature = "dice"))]
pub mod dice;

#[cfg(any(test, feature = "asserts"))]
pub mod asserts;

#[cfg(test)]
mod test_util;
