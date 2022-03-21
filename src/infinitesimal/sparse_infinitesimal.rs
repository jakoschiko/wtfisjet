use intmap::IntMap;

use crate::{Dim, Infinitesimal, Number};

/// Implementation for [`Infinitesimal`] that uses a sparse representation.
///
/// # Supported dimension count
///
/// Because of the limitations of the underlying data structure, the dimension count must fit
/// into `u64`. Any attempt to create an instance with greater dimension count will cause panic.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SparseInfinitesimal<N: Number> {
    cast_dim: u64, // Necessary because elems contains only non-zero elements
    zero: N,       // Necessary for functions that return `&N`
    elems: IntMap<N>,
}

impl<N: Number> Infinitesimal<N> for SparseInfinitesimal<N> {
    type DenseElems<'a> = SparseInfinitesimalDenseElems<'a, N> where N: 'a;

    type SparseElems<'a> = SparseInfinitesimalSparseElems<'a, N> where N: 'a;

    fn dim(&self) -> Dim {
        // Should never fail because the original `dim` value was of type `usize`.
        Dim(self.cast_dim as usize)
    }

    fn zeros(dim: Dim) -> Self {
        if let Ok(cast_dim) = dim.0.try_into() {
            Self {
                cast_dim,
                zero: N::zero(),
                elems: IntMap::new(),
            }
        } else {
            panic!(
                "SparseInfinitesimal doesn't support dimension count {0} because it doesn't fit into u64",
                dim.0
            )
        }
    }

    fn one(idx: usize, dim: Dim) -> Self {
        if idx >= dim.0 {
            panic!(
                "Index {idx} must not be equal to or greater than dimension count {}",
                dim.0
            )
        } else {
            // Should never fail because we already checked that `dim` fits into `u64` and
            // `idx` is smaller than `dim`.
            let cast_idx = idx as u64;
            let mut result = Self::zeros(dim);
            result.elems.insert(cast_idx, N::one());
            result
        }
    }

    fn from_dense<E: IntoIterator<Item = N>>(dense_elems: E) -> Self {
        let iter = dense_elems.into_iter();
        let mut elems = IntMap::with_capacity(iter.size_hint().0);
        let mut dim = Dim(0);

        for elem in iter {
            if !elem.is_zero() {
                // It doesn't matter that this could fail, we'll check `dim` later
                let cast_idx = dim.0 as u64;
                elems.insert(cast_idx, elem);
            }

            dim.0 += 1;
        }

        if let Ok(cast_dim) = dim.0.try_into() {
            Self {
                cast_dim,
                zero: N::zero(),
                elems,
            }
        } else {
            panic!(
                "SparseInfinitesimal doesn't support dimension count {} because it doesn't fit into u64",
                dim.0
            )
        }
    }

    fn from_sparse<E: IntoIterator<Item = (usize, N)>>(sparse_elems: E, dim: Dim) -> Self {
        let mut elems = IntMap::with_capacity(dim.0);

        let cast_dim = if let Ok(cast_dim) = dim.0.try_into() {
            cast_dim
        } else {
            panic!(
                "SparseInfinitesimal doesn't support dimension count {} because it doesn't fit into u64",
                dim.0
            )
        };

        for (idx, elem) in sparse_elems {
            let cast_idx = idx.try_into().ok().filter(|i| *i < cast_dim);

            if let Some(cast_idx) = cast_idx {
                if !elem.is_zero() {
                    // Unfortunately `insert` doesn't override the previous value, so we need
                    // to remove it first.
                    elems.remove(cast_idx);
                    elems.insert(cast_idx, elem);
                } else {
                    // In case of index collisions, we prefer the later value and the later value
                    // is zero, so we must take care that no value is stored for the index
                    elems.remove(cast_idx);
                }
            } else {
                panic!(
                    "Index {idx} must not be equal to or greater than dimension count {}",
                    dim.0
                )
            }
        }

        Self {
            cast_dim,
            zero: N::zero(),
            elems,
        }
    }

    fn is_sparse() -> bool {
        true
    }

    fn dense_elems(&self) -> Self::DenseElems<'_> {
        SparseInfinitesimalDenseElems {
            cast_idx: 0,
            cast_dim: self.cast_dim,
            zero: &self.zero,
            elems: &self.elems,
        }
    }

    fn sparse_elems(&self) -> Self::SparseElems<'_> {
        SparseInfinitesimalSparseElems {
            elems: self.elems.iter(),
        }
    }

    fn elem(&self, idx: usize) -> &N {
        let dim = self.dim();
        if idx >= dim.0 {
            panic!(
                "Index {idx} must not be equal to or greater than dimension count {}",
                dim.0
            )
        } else {
            // Should never fail because we already checked that `dim` fits into `u64` and
            // `idx` is smaller than `dim`.
            let cast_idx = idx as u64;
            self.elems.get(cast_idx).unwrap_or(&self.zero)
        }
    }

    fn add(mut self, rhs: Self) -> Self {
        let left_dim = self.cast_dim;
        let right_dim = rhs.cast_dim;

        if left_dim != right_dim {
            panic!(
                "Cannot add infinitesimal parts with different dimension counts {} and {}",
                left_dim, right_dim
            )
        } else {
            for (idx, right_elem) in rhs.elems.into_iter() {
                // Unfortunately, there is no entry API for `IntMap`.
                // So we need to check the key twice.
                if self.elems.contains_key(idx) {
                    let left_elem = self.elems.get_mut(idx).unwrap();
                    *left_elem += right_elem;

                    if left_elem.is_zero() {
                        self.elems.remove(idx);
                    }
                } else {
                    self.elems.insert(idx, right_elem);
                }
            }
            self
        }
    }

    fn sub(mut self, rhs: Self) -> Self {
        let left_dim = self.cast_dim;
        let right_dim = rhs.cast_dim;

        if left_dim != right_dim {
            panic!("Cannot subtract infinitesimal parts with different dimension counts {left_dim} and {right_dim}")
        } else {
            for (idx, right_elem) in rhs.elems.into_iter() {
                // Unfortunately, there is no entry API for `IntMap`.
                // So we need to check the key twice.
                if self.elems.contains_key(idx) {
                    let left_elem = self.elems.get_mut(idx).unwrap();
                    *left_elem -= right_elem;

                    if left_elem.is_zero() {
                        self.elems.remove(idx);
                    }
                } else {
                    self.elems.insert(idx, right_elem.neg());
                }
            }
            self
        }
    }

    fn neg(mut self) -> Self {
        for elem in self.elems.values_mut() {
            elem.neg_in_place();

            // We expect that `neg_in_place` doesn't turn a non-zero value into a zero value,
            // hence we don't need to check for zero here.
        }
        self
    }

    fn scale(mut self, factor: N) -> Self {
        if factor.is_zero() {
            // We expect that `MulAssign::mul_assign` applied to zero will produce only
            // zeros, therefore we implement a shortcut.
            self.elems.clear()
        } else {
            for elem in self.elems.values_mut() {
                *elem *= factor.clone()
            }

            // Unfortunately, the retain operator on `IntMap` doesn't emit mutable references,
            // so we need to iterate twice.
            self.elems.retain(|_, elem| !elem.is_zero())
        }
        self
    }
}

pub struct SparseInfinitesimalDenseElems<'a, N: 'a + Number> {
    cast_idx: u64,
    cast_dim: u64,
    zero: &'a N,
    elems: &'a IntMap<N>,
}

impl<'a, N: 'a + Number> Iterator for SparseInfinitesimalDenseElems<'a, N> {
    type Item = &'a N;

    fn next(&mut self) -> Option<Self::Item> {
        if self.cast_idx >= self.cast_dim {
            None
        } else {
            let result = self.elems.get(self.cast_idx).unwrap_or(self.zero);
            self.cast_idx += 1;
            Some(result)
        }
    }
}

pub struct SparseInfinitesimalSparseElems<'a, N: 'a + Number> {
    elems: intmap::Iter<'a, u64, N>,
}

impl<'a, N: Number> Iterator for SparseInfinitesimalSparseElems<'a, N> {
    type Item = (usize, &'a N);

    fn next(&mut self) -> Option<Self::Item> {
        self.elems.next().map(|(&idx, elem)| {
            // Should never fail because the original `idx` values were of type `usize`
            let cast_idx = idx as usize;
            (cast_idx, elem)
        })
    }
}

#[cfg(test)]
mod tests {
    use dicetest::prelude::*;

    use crate::SparseInfinitesimal;

    #[test]
    fn valid_infinitesimal_impl() {
        let dim_die = crate::dice::dim();
        let numbers_die = crate::dice::big_rational_non_zero_number();

        crate::asserts::assert_valid_infinitesimal_impl::<_, SparseInfinitesimal<_>, _, _>(
            Dicetest::repeatedly(),
            dim_die,
            numbers_die,
        );
    }
}
