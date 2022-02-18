use crate::Number;

/// The n-dimensional infinitesimal part of a [`Jet`] that represents the derivatives.
///
/// We use a trait because it allows us to choose an optimal implementation for each use case.
/// For example:
/// - If the infinitesimal part contains many zeros, we can reduce the memory usage and
/// calculation time by using an implementation based on sparse vectors.
/// - If we know the dimension at compile time, we can avoid heap allocations by using an
/// implementation based on arrays.
/// - If we want to reuse a function based on [`Jet`] but we are not interested in the derivatives,
/// we can prevent all calculations on the infinitesimal part by using a dummy implementation
/// with zero dimension.
///
/// [`Jet`]: crate::Jet
pub trait Infinitesimal<N: Number>: Clone + PartialEq {}
