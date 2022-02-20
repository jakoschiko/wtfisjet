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
pub trait Infinitesimal<N: Number>: Clone + PartialEq {
    /// Return type of [`Infinitesimal::dense_elems`].
    type DenseIterator: Iterator<Item = N>;

    /// Return type of [`Infinitesimal::sparse_elems`].
    type SparseIterator: Iterator<Item = (usize, N)>;

    /// Returns the dimension count of the infinitesimal part.
    fn dim(&self) -> usize;

    /// Return an instance with a one for the given dimension and zero for all other dimensions.
    ///
    /// If we want to differentiate with respect to a variable, we can use this function to
    /// initialize the infinitesimal part.
    ///
    /// # Panic
    ///
    /// This function panics if the given index is equal to or greater than the dimension count.
    fn one(idx: usize, dim: usize) -> Self;

    /// Returns an instance with zero for all dimensions.
    ///
    /// If want to create a constant jet without dependencies to any variables, we can use this
    /// function to initialize its infinitesimal part.
    fn zeros(dim: usize) -> Self;

    /// Returns whether the implementation uses a sparse representation that omits the zero
    /// elements.
    ///
    /// You can use this information for optimizations, e.g. store the non-zero elements
    /// inside a sparse matrix instead of dense matrix.
    fn is_sparse() -> bool;

    /// Returns the elements of the infinitesimal part for all dimensions.
    fn dense_elems(&self) -> Self::DenseIterator;

    /// Returns the non-zero elements of the infinitesimal part along with their dimensions.
    ///
    /// Each dimension must be emitted at most once. The order is undefined.
    fn sparse_elems(&self) -> Self::DenseIterator;

    /// Returns the element for the given dimension.
    ///
    /// # Panic
    ///
    /// This function panics if the given index is equal to or greater than the dimension count.
    fn elem(&self, idx: usize) -> N;

    /// Adds both infinitesimal parts element wise.
    ///
    /// # Panic
    ///
    /// This function panics if the infinitesimal parts have different dimension counts.
    fn add(self, lhs: Self) -> Self;

    /// Subtracts both infinitesimal parts element wise.
    ///
    /// # Panic
    ///
    /// This function panics if the infinitesimal parts have different dimension counts.
    fn sub(self, lhs: Self) -> Self;

    /// Negates all elements of the infinitesimal part.
    fn neg(self) -> Self;

    /// Scales all elements of the infinitesimal part with the given factor.
    fn scale(self, factor: N) -> Self;
}
