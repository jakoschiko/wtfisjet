use std::collections::BTreeMap;
use std::panic::{RefUnwindSafe, UnwindSafe};
use std::{collections::BTreeSet, fmt::Debug};

use diceprop::{ops, Elem, Fun1, Fun2, Set, Vars};
use dicetest::hint_section;
use dicetest::prelude::*;

use crate::{Dim, Infinitesimal, Number};

/// Asserts that the given implementation of [`Infinitesimal`] fulfills all required properties.
pub fn assert_valid_infinitesimal_impl<N, I, DD, DN>(
    dicetest: Dicetest,
    dim_die: DD,
    number_die: DN,
) where
    N: Number + PartialEq,
    I: Infinitesimal<N> + Debug,
    DD: Die<Dim> + RefUnwindSafe + UnwindSafe,
    DN: Die<N> + RefUnwindSafe + UnwindSafe,
{
    dicetest.run(move |mut fate| {
        let dim = fate.roll(&dim_die);
        hint!("dim = {:?}", dim);

        let valid_idx_die = if dim.0 == 0 { None } else { Some(dice::uni_usize(dim.indices())) };
        let infinitesimal_die = any_infinitesimal::<_, I, _>(dim, &number_die);
        let infinitesimal = fate.roll(&infinitesimal_die);

        // We assume that `I::dim` and `I::elem` are correctly implemented. The test expectations
        // of all other functions are based on these two functions.

        {
            hint_section!("Is `I::zeros` lawful?");

            let infinitesimal = I::zeros(dim);
            hint!("`I::zeros()` = {infinitesimal:?}");

            let infinitesimal_dim = infinitesimal.dim();
            assert_eq!(infinitesimal_dim, dim, "Result of `I::zeros` has wrong dimension {}", infinitesimal_dim.0);

            for idx in dim.indices() {
                let elem = infinitesimal.elem(idx);
                assert!(elem.is_zero(), "The element of `I::zeros` with index {idx} must be zero, but is {elem:?}");
            }
        }

        {
            hint_section!("Is `I::one` lawful?");

            if let Some(ref valid_idx_die) = valid_idx_die {
                let one_idx = fate.roll(valid_idx_die);

                let infinitesimal = I::one(one_idx, dim);
                hint!("`I::one({one_idx}, {})` = {infinitesimal:?}", dim.0);

                let infinitesimal_dim = infinitesimal.dim();
                assert_eq!(infinitesimal_dim, dim, "Result of `I::one` has wrong dimension {}", infinitesimal_dim.0);

                for idx in dim.indices() {
                    let elem = infinitesimal.elem(idx);

                    if idx == one_idx {
                        assert!(elem.is_one(), "The element of `I::one` with index {idx} must be one, but is {elem:?}");
                    } else {
                        assert!(elem.is_zero(), "The element of `I::one` with index {idx} must be zero, but is {elem:?}");
                    }
                }
            }
        }

        {
            hint_section!("Is `I::from_dense` lawful?");

            let numbers = fate.roll(dice::vec(&number_die, dim.0));

            let infinitesimal = I::from_dense(numbers.clone());
            hint!("`I::from_dense({numbers:?})` = {infinitesimal:?}");

            let infinitesimal_dim = infinitesimal.dim();
            assert_eq!(infinitesimal_dim, dim, "Result of `I::from_dense` has wrong dimension {}", infinitesimal_dim.0);

            for idx in dim.indices() {
                let elem = infinitesimal.elem(idx);
                let expected = &numbers[idx];

                assert_eq!(elem, expected, "The element of `I::from_dense` with index {idx} must be {expected:?}, but is {elem:?}");
            }
        }

        {
            hint_section!("Is `I::from_sparse` lawful?");

            if let Some(ref valid_idx_die) = valid_idx_die {
                let number_with_idx_die = dice::zip().two(
                    valid_idx_die,
                    // We want to generate zero with higher probability to reveal
                    // sparse related bugs.
                    dice::one_of_die().two(&number_die, dice::just(N::zero())),
                );
                let numbers = fate.roll(dice::vec(number_with_idx_die, ..));
                // In case of index collisions, the later number will be taken
                let mut spare_numbers_map = numbers.clone().into_iter().collect::<BTreeMap<_,_>>();
                spare_numbers_map.retain(|_, n| !n.is_zero());

                let infinitesimal = I::from_sparse(numbers.clone(), dim);
                hint!("`I::from_sparse({numbers:?}, dim)` = {infinitesimal:?}");

                let infinitesimal_dim = infinitesimal.dim();
                assert_eq!(infinitesimal_dim, dim, "Result of `I::from_sparse` has wrong dimension {}", infinitesimal_dim.0);

                for idx in dim.indices() {
                    let elem = infinitesimal.elem(idx);

                    if let Some(expected) = spare_numbers_map.get(&idx) {
                        assert_eq!(elem, expected, "The element of `I::from_sparse` with index {idx} must be {expected:?}, but is {elem:?}");
                    } else {
                        assert!(elem.is_zero(), "The element of `I::from_sparse` with index {idx} must be zero, but is {elem:?}");
                    }
                }
            }
        }

        {
            hint_section!("Is `I::dense_elems` lawful?");

            let dense_elems = infinitesimal.dense_elems().cloned().collect::<Vec<_>>();
            hint!("`I::dense_elems({infinitesimal:?})` = {dense_elems:?}");

            let dense_elems_len = dense_elems.len();
            assert_eq!(dense_elems_len, dim.0, "Result of `I::dense_elems` has wrong length {dense_elems_len}");

            for idx in dim.indices() {
                let elem = &dense_elems[idx];
                let expected = infinitesimal.elem(idx);

                assert_eq!(elem, expected, "The element of `I::dense_elems` with index {idx} must be {expected:?}, but is {elem:?}");
            }
        }

        {
            hint_section!("Is `I::sparse_elems` lawful?");

            let sparse_elems = infinitesimal.sparse_elems().collect::<Vec<_>>();
            hint!("`I::sparse_elems({infinitesimal:?})` = {sparse_elems:?}");

            let mut seen_indices = BTreeSet::new();

            for (idx, elem) in sparse_elems {
                assert!(idx < dim.0, "`I::sparse_elems` must emit only indices smaller than the dimension count, but index {idx} was emitted");

                assert!(!seen_indices.contains(&idx), "`I::sparse_elems` must emit each index at most once, but the index {idx} was emitted multiple times");
                seen_indices.insert(idx);

                assert!(!elem.is_zero(), "`I::sparse_elems` must not emit zero, but zero was emitted for index {idx}");

                let expected = infinitesimal.elem(idx);

                assert_eq!(elem, expected, "The element of `I::sparse_elems` with index {idx} must be {expected:?}, but is {elem:?}");
            }

            let mut expected_indices = BTreeSet::new();

            for idx in dim.indices() {
                if !infinitesimal.elem(idx).is_zero() {
                    expected_indices.insert(idx);
                }
            }

            assert_eq!(seen_indices, expected_indices, "The indices emitted by `I::sparse_elems` must be {expected_indices:?}, but is {seen_indices:?}");
        }

        {
            hint_section!("Are `I::add`, `I::sub`, `I::neg`, `I::scale` lawful?");

            let i_set = Set::new("I", &infinitesimal_die);
            let Vars { set: i_set_name, elems: [x, y, z] } = fate.roll(i_set.vars(["x", "y", "z"]));
            let i_vars1 = Vars::new(i_set_name, [x.clone()]);
            let i_vars2 = Vars::new(i_set_name, [x.clone(), y.clone()]);
            let i_vars3 = Vars::new(i_set_name, [x, y, z]);
            let i_zero = Elem::new("zero", I::zeros(dim));
            let i_add = Fun2::new("I::add", I::add);
            let i_sub = Fun2::new("I::sub", I::sub);
            let i_neg = Fun1::new("I::neg", I::neg);

            let n_set = Set::new("N", &number_die);
            let n_vars2 = fate.roll(n_set.vars(["a", "b"]));

            diceprop::props::algebra::abelian_group(i_vars3, i_add.as_ref(), i_neg, i_zero);
            diceprop::props::binop::inverse(i_vars2, i_add, i_sub);
            prop_left_distributive(i_vars1, n_vars2);
        }
    })
}

// Generator for `Infinitesimal` that is based on any available function that returns
// an instance.
//
// This is useful because in combination with the test above this might reveal
// that one of these functions returns an invalid instance, e.g. an instance
// with sparse representation that contains a zero.
fn any_infinitesimal<N: Number, I: Infinitesimal<N>, NDI: Die<N>>(
    dim: Dim,
    infinitesimal_number_die: NDI,
) -> impl Die<I> {
    dice::from_fn(move |mut fate| {
        // This already creates instances via `I::zeroes`, `I::one`, `I::from_dense` and
        // `I::from_sparse`
        let infinitesimal_constructor_die =
            crate::dice::infinitesimal(dim, &infinitesimal_number_die);

        let infinitesimal_add_die = dice::zip()
            .two(
                &infinitesimal_constructor_die,
                &infinitesimal_constructor_die,
            )
            .map(|(l, r)| I::add(l, r));
        let infinitesimal_sub_die = dice::zip()
            .two(
                &infinitesimal_constructor_die,
                &infinitesimal_constructor_die,
            )
            .map(|(l, r)| I::sub(l, r));
        let infinitesimal_neg_die = (&infinitesimal_constructor_die).map(I::neg);
        let infinitesimal_scale_die = dice::zip()
            .two(&infinitesimal_constructor_die, &infinitesimal_number_die)
            .map(|(i, f)| I::scale(i, f));

        let infinitesimal_die = dice::one_of_die().five(
            &infinitesimal_constructor_die,
            &infinitesimal_add_die,
            &infinitesimal_sub_die,
            &infinitesimal_neg_die,
            &infinitesimal_scale_die,
        );

        fate.roll(infinitesimal_die)
    })
}

pub fn prop_left_distributive<N, I>(i_var: Vars<I, 1>, n_vars: Vars<N, 2>)
where
    N: Number + PartialEq,
    I: Infinitesimal<N> + Debug,
{
    let i_add = Fun2::new("I::add", I::add);
    let n_add = Fun2::new("N::add", N::add);
    let i_scale = Fun2::new("I::scale", I::scale);

    hint_section!(
        "Is `{}` left distributive over `{}`?",
        i_scale.name,
        i_add.name,
    );

    let [x] = i_var.eval();
    let [a, b] = n_vars.eval();

    ops::assert(ops::eq(
        i_scale
            .eval(x.clone(), n_add.eval(a.clone(), b.clone()))
            .as_ref(),
        i_add
            .eval(
                i_scale.eval(x.clone(), a.clone()),
                i_scale.eval(x.clone(), b.clone()),
            )
            .as_ref(),
    ));
}
