use std::marker::PhantomData;

use crate::{Infinitesimal, Number};

/// Implementation for [`Infinitesimal`] with zero dimensions.
///
/// Useful if you want to use a jet based function but you are not interested
/// in the derivatives. Because it's a zero size type with only no-ops, all
/// overhead caused the infinitesimal part of the jet will be eliminated by
/// the compiler.
///
/// # Supported dimension count
///
/// This implementation supports only the dimension count zero. Any attempt to
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

    fn zeros(dim: usize) -> Self {
        if dim == 0 {
            NoInfinitesimal
        } else {
            panic!("NoInfinitesimal doesn't support dimension count {dim}")
        }
    }

    fn one(idx: usize, dim: usize) -> Self {
        if dim == 0 {
            panic!("NoInfinitesimal doesn't support dimension with index {idx}")
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
