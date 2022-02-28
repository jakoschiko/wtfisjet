//! Generators for random test data.

use dicetest::prelude::*;

use crate::{Dim, Infinitesimal, Jet, Number};

#[cfg(any(test, feature = "big-rational-number"))]
/// Generates an arbitrary [`BigRational`].
pub fn big_rational_number() -> impl Die<num_rational::BigRational> {
    use num_traits::{One, Zero};

    let special_number_die = dice::one_of().three(
        num_rational::BigRational::zero(),
        num_rational::BigRational::one(),
        -num_rational::BigRational::one(),
    );

    let number_die = dice::from_fn(|mut fate| {
        let numerator = fate.roll(dice::u32(..)).into();
        let denominator = fate.roll(dice::i32(1..)).into();
        num_rational::BigRational::new(numerator, denominator)
    });

    dice::weighted_one_of_die().two((1, special_number_die), (7, number_die))
}

#[cfg(any(test, feature = "big-rational-number"))]
/// Generates an arbitrary non-zero [`BigRational`].
pub fn big_rational_non_zero_number() -> impl Die<num_rational::BigRational> {
    use num_traits::One;

    let special_number_die = dice::one_of().two(
        num_rational::BigRational::one(),
        -num_rational::BigRational::one(),
    );

    let number_die = dice::from_fn(|mut fate| {
        let numerator = fate.roll(dice::u32(1..)).into();
        let denominator = fate.roll(dice::i32(1..)).into();
        num_rational::BigRational::new(numerator, denominator)
    });

    dice::weighted_one_of_die().two((1, special_number_die), (7, number_die))
}

/// Generates an arbitrary [`Dim`] that is limited by [`dicetest::Limit`].
pub fn dim() -> impl Die<Dim> {
    dice::length(..).map(Dim)
}

/// Generates an [`Infinitesimal`] using [`Infinitesimal::one`].
///
/// # Panic
///
/// This function panics if the given dimension count is zero.
pub fn infinitesimal_one<N: Number, I: Infinitesimal<N>>(dim: Dim) -> impl Die<I> {
    assert!(
        dim.0 != 0,
        "Generator infinitesimal_one must not be used with dimension count 0"
    );

    dice::from_fn(move |mut fate| {
        let idx = fate.roll(dice::uni_usize(0..dim.0));
        I::one(idx, dim)
    })
}

/// Generates an [`Infinitesimal`] with arbitrary elements.
pub fn infinitesimal<N: Number, I: Infinitesimal<N>, NDI: Die<N>>(
    dim: Dim,
    infinitesimal_number_die: NDI,
) -> impl Die<I> {
    dice::from_fn(move |mut fate| {
        let infinitesimal_zeros_die = dice::from_fn(|_| I::zeros(dim));

        let safe_infinitesimal_one_die = dice::from_fn(|mut fate| {
            if dim.0 == 0 {
                I::zeros(Dim(0))
            } else {
                fate.roll(infinitesimal_one(dim))
            }
        });

        let infinitesimal_dense_die = dice::from_fn(|mut fate| {
            let elems = fate.roll(dice::vec(&infinitesimal_number_die, dim.0));
            I::from_dense(elems)
        });

        let infinitesimal_sparse_die = dice::from_fn(|mut fate| {
            if dim.0 == 0 {
                I::zeros(Dim(0))
            } else {
                let number_with_index_die =
                    dice::zip().two(dice::uni_usize(0..dim.0), &infinitesimal_number_die);
                let elems = fate.roll(dice::vec(number_with_index_die, 0..=dim.0));
                I::from_sparse(elems, dim)
            }
        });

        // Compose different generators because each of them has a different distribution
        let infinitesimal_die = dice::one_of_die().three(
            dice::one_of_die().two(&infinitesimal_zeros_die, &safe_infinitesimal_one_die),
            &infinitesimal_dense_die,
            &infinitesimal_sparse_die,
        );

        fate.roll(infinitesimal_die)
    })
}

/// Generates a [`Jet`] using [`Jet::constant`].
pub fn jet_constant<N: Number, I: Infinitesimal<N>, NDR: Die<N>>(
    dim: Dim,
    real_number_die: NDR,
) -> impl Die<Jet<N, I>> {
    dice::from_fn(move |mut fate| {
        let real = fate.roll(&real_number_die);
        Jet::constant(real, dim)
    })
}

/// Generates a [`Jet`] using [`Jet::variable`].
///
/// # Panic
///
/// This function panics if the given dimension count is zero.
pub fn jet_variable<N: Number, I: Infinitesimal<N>, NDR: Die<N>>(
    dim: Dim,
    real_number_die: NDR,
) -> impl Die<Jet<N, I>> {
    assert!(
        dim.0 != 0,
        "Generator jet_variable must not be used with dimension count 0"
    );

    dice::from_fn(move |mut fate| {
        let real = fate.roll(&real_number_die);
        let idx = fate.roll(dice::uni_usize(0..dim.0));
        Jet::variable(real, idx, dim)
    })
}

/// Generates a [`Jet`] with arbitrary values for the real part and the infinitesimal part.
pub fn jet<N: Number, I: Infinitesimal<N>, NDR: Die<N>, NDI: Die<N>>(
    dim: Dim,
    real_number_die: NDR,
    infinitesimal_number_die: NDI,
) -> impl Die<Jet<N, I>> {
    dice::from_fn(move |mut fate| {
        let any_jet_die = dice::zip()
            .two(
                &real_number_die,
                infinitesimal(dim, &infinitesimal_number_die),
            )
            .map_once(|(real, infinitesimal)| Jet::new(real, infinitesimal));

        let safe_variable_die = dice::from_fn(|mut fate| {
            if dim.0 == 0 {
                fate.roll(jet_constant(dim, &real_number_die))
            } else {
                fate.roll(jet_variable(dim, &real_number_die))
            }
        });

        // Compose different generators because each of them has a different distribution
        let jet_die = dice::one_of_die_once().three(
            jet_constant(dim, &real_number_die),
            safe_variable_die,
            any_jet_die,
        );
        fate.roll(jet_die)
    })
}
