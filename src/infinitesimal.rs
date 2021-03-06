use std::fmt::Debug;

use crate::{Dim, Number};

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
pub trait Infinitesimal<N: Number>: Debug + Clone + PartialEq {
    /// Return type of [`Infinitesimal::dense_elems`].
    type DenseElems<'a>: Iterator<Item = &'a N>
    where
        Self: 'a,
        N: 'a;

    /// Return type of [`Infinitesimal::sparse_elems`].
    type SparseElems<'a>: Iterator<Item = (usize, &'a N)>
    where
        Self: 'a,
        N: 'a;

    /// Returns the dimension count of the infinitesimal part.
    fn dim(&self) -> Dim;

    /// Returns an instance with zero for all dimensions.
    ///
    /// If want to create a constant jet without dependencies to any variables, we can use this
    /// function to initialize its infinitesimal part.
    ///
    /// # Panic
    ///
    /// This function panics if the implementation does not support the dimension count.
    fn zeros(dim: Dim) -> Self;

    /// Return an instance with a one for the given dimension and zero for all other dimensions.
    ///
    /// If we want to differentiate with respect to a variable, we can use this function to
    /// initialize the infinitesimal part.
    ///
    /// # Panic
    ///
    /// This function panics if
    /// - the given index is equal to or greater than the dimension count or
    /// - the implementation does not support the dimension count.
    fn one(idx: usize, dim: Dim) -> Self;

    /// Returns an instance that contains the elements emitted by the given iterator.
    ///
    /// The number of elements emitted by the iterator is the dimension count of the
    /// returned instance.
    ///
    /// # Panic
    ///
    /// This function panics if the implementation does not support the dimension count.
    fn from_dense<E: IntoIterator<Item = N>>(dense_elems: E) -> Self;

    /// Returns an instance that contains the elements emitted by the given iterator
    /// and zero for all other dimensions.
    ///
    /// Each element is emitted along with its dimension index. The indices are allowed
    /// to be emitted in arbitrary order. If an index is emitted multiple times, the last
    /// element will be used.
    ///
    /// # Panic
    ///
    /// This function panics if
    /// - the iterator emits an element with an index equal to or greater than the
    /// dimension count or
    /// - the implementation does no support the dimension count.
    fn from_sparse<E: IntoIterator<Item = (usize, N)>>(sparse_elems: E, dim: Dim) -> Self;

    /// Returns whether the implementation uses a sparse representation that omits the zero
    /// elements.
    ///
    /// You can use this information for optimizations, e.g. store the non-zero elements
    /// inside a sparse matrix instead of dense matrix.
    fn is_sparse() -> bool;

    /// Returns the elements of the infinitesimal part for all dimensions.
    fn dense_elems(&self) -> Self::DenseElems<'_>;

    /// Returns the non-zero elements of the infinitesimal part along with their dimensions.
    ///
    /// Each dimension must be emitted at most once. The order is undefined.
    fn sparse_elems(&self) -> Self::SparseElems<'_>;

    /// Returns the element for the given dimension.
    ///
    /// # Panic
    ///
    /// This function panics if the given index is equal to or greater than the dimension count.
    fn elem(&self, idx: usize) -> &N;

    /// Adds both infinitesimal parts element wise.
    ///
    /// # Panic
    ///
    /// This function panics if the infinitesimal parts have different dimension counts.
    fn add(self, rhs: Self) -> Self;

    /// Subtracts both infinitesimal parts element wise.
    ///
    /// # Panic
    ///
    /// This function panics if the infinitesimal parts have different dimension counts.
    fn sub(self, rhs: Self) -> Self;

    /// Negates all elements of the infinitesimal part.
    fn neg(self) -> Self;

    /// Scales all elements of the infinitesimal part with the given factor.
    fn scale(self, factor: N) -> Self;
}

mod no_infinitesimal;
pub use no_infinitesimal::NoInfinitesimal;

mod dense_infinitesimal;
pub use dense_infinitesimal::DenseInfinitesimal;

#[cfg(any(test, feature = "sparse-infinitesimal"))]
mod sparse_infinitesimal;
#[cfg(any(test, feature = "sparse-infinitesimal"))]
pub use sparse_infinitesimal::SparseInfinitesimal;

#[cfg(any(test, feature = "const-infinitesimal"))]
mod const_infinitesimal;
#[cfg(any(test, feature = "const-infinitesimal"))]
pub use const_infinitesimal::ConstInfinitesimal;
