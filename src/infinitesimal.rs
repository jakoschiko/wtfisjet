use std::marker::PhantomData;

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
    fn dim(&self) -> usize;

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
    fn one(idx: usize, dim: usize) -> Self;

    /// Returns an instance with zero for all dimensions.
    ///
    /// If want to create a constant jet without dependencies to any variables, we can use this
    /// function to initialize its infinitesimal part.
    ///
    /// # Panic
    ///
    /// This function panics if the implementation does no support the dimension count.
    fn zeros(dim: usize) -> Self;

    /// Returns an instance that contains the elements emitted by the given iterator.
    ///
    /// The number of elements emitted by the iterator is the dimension count of the
    /// returned instance.
    ///
    /// # Panic
    ///
    /// This function panics if the implementation does no support the dimension count.
    fn from_dense<I: Iterator<Item = N>>(dense_elems: I) -> Self;

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
    fn from_sparse<I: Iterator<Item = (usize, N)>>(sparse_elems: I, dim: usize) -> Self;

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

/// Implementation for [`Infinitesimal`] with zero dimensions.
///
/// Useful if you want to use a jet based function but you are not interested
/// in the derivatives. Because it's a zero size type with only no-ops, all
/// overhead caused the infinitesimal part of the jet will be eliminated by
/// the compiler.
///
/// # Supported dimension count
///
/// This implementation supports only the dimension count zero. All attempts to
/// create an instance with greater dimension count will cause panic.
#[derive(Debug, Clone, PartialEq)]
pub struct NoInfinitesimal;

impl<N: Number> Infinitesimal<N> for NoInfinitesimal {
    type DenseElems<'a>
    where
        N: 'a,
    = NoInfinitesimalDenseElems<'a, N>;

    type SparseElems<'a>
    where
        N: 'a,
    = NoInfinitesimalSparseElems<'a, N>;

    fn dim(&self) -> usize {
        0
    }

    fn one(idx: usize, dim: usize) -> Self {
        if dim == 0 {
            panic!("NoInfinitesimal doesn't support dimension with index {idx}")
        } else {
            panic!("NoInfinitesimal doesn't support dimension count {dim}")
        }
    }

    fn zeros(dim: usize) -> Self {
        if dim == 0 {
            NoInfinitesimal
        } else {
            panic!("NoInfinitesimal doesn't support dimension count {dim}")
        }
    }

    fn from_dense<I: Iterator<Item = N>>(dense_elems: I) -> Self {
        let dim = dense_elems.count();
        if dim == 0 {
            NoInfinitesimal
        } else {
            panic!("NoInfinitesimal doesn't support dimension count {dim}")
        }
    }

    fn from_sparse<I: Iterator<Item = (usize, N)>>(mut sparse_elems: I, dim: usize) -> Self {
        if dim == 0 {
            if let Some((idx, _)) = sparse_elems.next() {
                panic!("NoInfinitesimal doesn't support dimension with index {idx}")
            } else {
                NoInfinitesimal
            }
        } else {
            panic!("NoInfinitesimal doesn't support dimension count {dim}")
        }
    }

    fn is_sparse() -> bool {
        false
    }

    fn dense_elems(&self) -> Self::DenseElems<'_> {
        NoInfinitesimalDenseElems { _n: PhantomData }
    }

    fn sparse_elems(&self) -> Self::SparseElems<'_> {
        NoInfinitesimalSparseElems { _n: PhantomData }
    }

    fn elem(&self, idx: usize) -> &N {
        panic!("NoInfinitesimal doesn't support dimension with index {idx}")
    }

    fn add(self, _rhs: Self) -> Self {
        NoInfinitesimal
    }

    fn sub(self, _rhs: Self) -> Self {
        NoInfinitesimal
    }

    fn neg(self) -> Self {
        NoInfinitesimal
    }

    fn scale(self, _factor: N) -> Self {
        NoInfinitesimal
    }
}

pub struct NoInfinitesimalDenseElems<'a, N: 'a + Number> {
    _n: PhantomData<&'a N>,
}

impl<'a, N: 'a + Number> Iterator for NoInfinitesimalDenseElems<'a, N> {
    type Item = &'a N;

    fn next(&mut self) -> Option<Self::Item> {
        None
    }
}

pub struct NoInfinitesimalSparseElems<'a, N: 'a + Number> {
    _n: PhantomData<&'a N>,
}

impl<'a, N: Number> Iterator for NoInfinitesimalSparseElems<'a, N> {
    type Item = (usize, &'a N);

    fn next(&mut self) -> Option<Self::Item> {
        None
    }
}

/// Implementation for [`Infinitesimal`] that uses a dense representation.
///
/// # Supported dimension count
///
/// There are no constraints regarding the dimension count.
#[derive(Debug, Clone, PartialEq)]
pub struct DenseInfinitesimal<N: Number> {
    elems: Vec<N>,
}

impl<N: Number> Infinitesimal<N> for DenseInfinitesimal<N> {
    type DenseElems<'a>
    where
        N: 'a,
    = DenseInfinitesimalDenseElems<'a, N>;

    type SparseElems<'a>
    where
        N: 'a,
    = DenseInfinitesimalSparseElems<'a, N>;

    fn dim(&self) -> usize {
        self.elems.len()
    }

    fn one(idx: usize, dim: usize) -> Self {
        if idx >= dim {
            panic!("Index {idx} must not be equal to or greater than dimension count {dim}")
        } else {
            let mut result = Self::zeros(dim);
            result.elems[idx] = N::one();
            result
        }
    }

    fn zeros(dim: usize) -> Self {
        let elems = std::iter::repeat_with(N::zero).take(dim).collect();
        Self { elems }
    }

    fn from_dense<I: Iterator<Item = N>>(dense_elems: I) -> Self {
        let elems = dense_elems.collect();
        Self { elems }
    }

    fn from_sparse<I: Iterator<Item = (usize, N)>>(sparse_elems: I, dim: usize) -> Self {
        let mut result = Self::zeros(dim);

        for (idx, elem) in sparse_elems {
            if idx >= dim {
                panic!("Index {idx} must not be equal to or greater than dimension count {dim}")
            } else {
                result.elems[idx] = elem;
            }
        }

        result
    }

    fn is_sparse() -> bool {
        false
    }

    fn dense_elems(&self) -> Self::DenseElems<'_> {
        DenseInfinitesimalDenseElems {
            elems: self.elems.iter(),
        }
    }

    fn sparse_elems(&self) -> Self::SparseElems<'_> {
        DenseInfinitesimalSparseElems {
            idx: 0,
            elems: self.elems.iter(),
        }
    }

    fn elem(&self, idx: usize) -> &N {
        let dim = self.dim();
        if idx >= dim {
            panic!("Index {idx} must not be equal to or greater than dimension count {dim}")
        } else {
            &self.elems[idx]
        }
    }

    fn add(mut self, rhs: Self) -> Self {
        let left_dim = self.dim();
        let right_dim = rhs.dim();

        if left_dim != right_dim {
            panic!("Cannot add infinitesimal parts with different dimension counts {left_dim} and {right_dim}")
        } else {
            for (left_elem, right_elem) in self.elems.iter_mut().zip(rhs.elems.into_iter()) {
                *left_elem += right_elem
            }
            self
        }
    }

    fn sub(mut self, rhs: Self) -> Self {
        let left_dim = self.dim();
        let right_dim = rhs.dim();

        if left_dim != right_dim {
            panic!("Cannot subtract infinitesimal parts with different dimension counts {left_dim} and {right_dim}")
        } else {
            for (left_elem, right_elem) in self.elems.iter_mut().zip(rhs.elems.into_iter()) {
                *left_elem -= right_elem
            }
            self
        }
    }

    fn neg(mut self) -> Self {
        for elem in self.elems.iter_mut() {
            // Unfortunately, there is no in-place version of `Neg`.
            // So we need to use a workaround.
            let temp = std::mem::replace(elem, N::zero());
            *elem = -temp;
        }
        self
    }

    fn scale(mut self, factor: N) -> Self {
        for elem in self.elems.iter_mut() {
            *elem *= factor.clone()
        }
        self
    }
}

pub struct DenseInfinitesimalDenseElems<'a, N: 'a + Number> {
    elems: std::slice::Iter<'a, N>,
}

impl<'a, N: 'a + Number> Iterator for DenseInfinitesimalDenseElems<'a, N> {
    type Item = &'a N;

    fn next(&mut self) -> Option<Self::Item> {
        self.elems.next()
    }
}

pub struct DenseInfinitesimalSparseElems<'a, N: 'a + Number> {
    idx: usize,
    elems: std::slice::Iter<'a, N>,
}

impl<'a, N: Number> Iterator for DenseInfinitesimalSparseElems<'a, N> {
    type Item = (usize, &'a N);

    fn next(&mut self) -> Option<Self::Item> {
        for elem in self.elems.by_ref() {
            let idx = self.idx;
            self.idx += 1;

            if !elem.is_zero() {
                return Some((idx, elem));
            }
        }
        None
    }
}
