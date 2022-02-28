use crate::{Dim, Infinitesimal, Number};

/// Implementation for [`Infinitesimal`] with const dimension count.
///
/// Because the implementation is based on arrays, it allows you to perform [`Jet`] calculations
/// without any heap allocation. Unfortunately this only works if the dimension count is known
/// at compile time.
///
/// [`Jet`]: crate::Jet
///
/// # Supported dimension count
///
/// This implementation only allows to create instances with a dimension count that matches the
/// const parameter `D`. Any attempt to create an instance with another dimension count
/// will cause panic.
#[derive(Debug, Clone, PartialEq)]
pub struct ConstInfinitesimal<N: Number, const D: usize> {
    elems: [N; D],
}

impl<N: Number, const D: usize> Infinitesimal<N> for ConstInfinitesimal<N, D> {
    type DenseElems<'a>
    where
        N: 'a,
    = ConstInfinitesimalDenseElems<'a, N>;

    type SparseElems<'a>
    where
        N: 'a,
    = ConstInfinitesimalSparseElems<'a, N>;

    fn dim(&self) -> Dim {
        Dim(D)
    }

    fn zeros(dim: Dim) -> Self {
        if dim.0 != D {
            panic!(
                "ConstInfinitesimal<_, {D}> doesn't support dimension count {}",
                dim.0
            )
        } else {
            let elems = array_init::array_init(|_| N::zero());
            ConstInfinitesimal { elems }
        }
    }

    fn one(idx: usize, dim: Dim) -> Self {
        if dim.0 != D {
            panic!(
                "ConstInfinitesimal<_, {D}> doesn't support dimension count {}",
                dim.0
            )
        } else if idx >= D {
            panic!("Index {idx} must not be equal to or greater than dimension count {D}")
        } else {
            let mut result = Self::zeros(Dim(D));
            result.elems[idx] = N::one();
            result
        }
    }

    fn from_dense<E: IntoIterator<Item = N>>(dense_elems: E) -> Self {
        if let Some(elems) = array_init::from_iter(dense_elems) {
            Self { elems }
        } else {
            panic!("ConstInfinitesimal<_, {D}> doesn't support dimension count greater than {D}")
        }
    }

    fn from_sparse<E: IntoIterator<Item = (usize, N)>>(sparse_elems: E, dim: Dim) -> Self {
        if dim.0 != D {
            panic!(
                "ConstInfinitesimal<_, {D}> doesn't support dimension count {}",
                dim.0
            )
        } else {
            let mut result = Self::zeros(Dim(D));

            for (idx, elem) in sparse_elems {
                if let Some(target) = result.elems.get_mut(idx) {
                    *target = elem
                } else {
                    panic!("Index {idx} must not be equal to or greater than dimension count {D}")
                }
            }

            result
        }
    }

    fn is_sparse() -> bool {
        false
    }

    fn dense_elems(&self) -> Self::DenseElems<'_> {
        ConstInfinitesimalDenseElems {
            elems: self.elems.iter(),
        }
    }

    fn sparse_elems(&self) -> Self::SparseElems<'_> {
        ConstInfinitesimalSparseElems {
            idx: 0,
            elems: self.elems.iter(),
        }
    }

    fn elem(&self, idx: usize) -> &N {
        if idx >= D {
            panic!("Index {idx} must not be equal to or greater than dimension count {D}")
        } else {
            &self.elems[idx]
        }
    }

    fn add(mut self, rhs: Self) -> Self {
        for (left_elem, right_elem) in self.elems.iter_mut().zip(rhs.elems.into_iter()) {
            *left_elem += right_elem
        }
        self
    }

    fn sub(mut self, rhs: Self) -> Self {
        for (left_elem, right_elem) in self.elems.iter_mut().zip(rhs.elems.into_iter()) {
            *left_elem -= right_elem
        }
        self
    }

    fn neg(mut self) -> Self {
        for elem in self.elems.iter_mut() {
            elem.neg_in_place();
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

pub struct ConstInfinitesimalDenseElems<'a, N: 'a + Number> {
    elems: std::slice::Iter<'a, N>,
}

impl<'a, N: 'a + Number> Iterator for ConstInfinitesimalDenseElems<'a, N> {
    type Item = &'a N;

    fn next(&mut self) -> Option<Self::Item> {
        self.elems.next()
    }
}

pub struct ConstInfinitesimalSparseElems<'a, N: 'a + Number> {
    idx: usize,
    elems: std::slice::Iter<'a, N>,
}

impl<'a, N: Number> Iterator for ConstInfinitesimalSparseElems<'a, N> {
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

    use crate::{ConstInfinitesimal, Dim};

    fn valid_infinitesimal_impl<const N: usize>(dicetest: Dicetest) {
        let dim_die = dice::just(Dim(N));
        let numbers_die = crate::dice::big_rational_non_zero_number();

        crate::asserts::assert_valid_infinitesimal_impl::<_, ConstInfinitesimal<_, N>, _, _>(
            dicetest,
            dim_die,
            numbers_die,
        );
    }

    #[test]
    fn valid_infinitesimal_impl_0() {
        valid_infinitesimal_impl::<0>(Dicetest::repeatedly());
    }

    #[test]
    fn valid_infinitesimal_impl_1() {
        valid_infinitesimal_impl::<1>(Dicetest::repeatedly());
    }

    #[test]
    fn valid_infinitesimal_impl_2() {
        valid_infinitesimal_impl::<2>(Dicetest::repeatedly());
    }

    #[test]
    fn valid_infinitesimal_impl_3() {
        valid_infinitesimal_impl::<3>(Dicetest::repeatedly());
    }

    #[test]
    fn valid_infinitesimal_impl_4() {
        valid_infinitesimal_impl::<4>(Dicetest::repeatedly());
    }

    #[test]
    fn valid_infinitesimal_impl_17() {
        valid_infinitesimal_impl::<17>(Dicetest::repeatedly());
    }

    #[test]
    fn valid_infinitesimal_impl_71() {
        valid_infinitesimal_impl::<71>(Dicetest::repeatedly());
    }
}
