use crate::{Infinitesimal, Number};

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

    fn zeros(dim: usize) -> Self {
        let elems = std::iter::repeat_with(N::zero).take(dim).collect();
        Self { elems }
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

    fn from_dense<E: IntoIterator<Item = N>>(dense_elems: E) -> Self {
        let elems = dense_elems.into_iter().collect();
        Self { elems }
    }

    fn from_sparse<E: IntoIterator<Item = (usize, N)>>(sparse_elems: E, dim: usize) -> Self {
        let mut result = Self::zeros(dim);

        for (idx, elem) in sparse_elems {
            if let Some(target) = result.elems.get_mut(idx) {
                *target = elem
            } else {
                panic!("Index {idx} must not be equal to or greater than dimension count {dim}")
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

#[cfg(test)]
mod tests {
    use dicetest::prelude::*;

    use crate::DenseInfinitesimal;

    #[test]
    fn valid_infinitesimal_impl() {
        let dim_die = dice::length(..);
        let numbers_die = crate::dice::big_rational_non_zero_number();

        crate::asserts::assert_valid_infinitesimal_impl::<_, DenseInfinitesimal<_>, _, _>(
            Dicetest::repeatedly(),
            dim_die,
            numbers_die,
        );
    }
}
