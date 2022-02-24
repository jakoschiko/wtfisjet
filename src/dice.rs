//! Generators for random test data.

use dicetest::prelude::*;

use crate::{Infinitesimal, Jet, Number};

#[cfg(any(test, feature = "big-rational-number"))]
/// Generates an arbitrary [`BigRational`].
pub fn big_rational_number() -> impl Die<num_rational::BigRational> {
    dice::from_fn(|mut fate| {
        let numerator = fate.roll(dice::u32(..)).into();
        let denominator = fate.roll(dice::i32(1..)).into();
        num_rational::BigRational::new(numerator, denominator)
    })
}

#[cfg(any(test, feature = "big-rational-number"))]
/// Generates an arbitrary non-zero [`BigRational`].
pub fn big_rational_non_zero_number() -> impl Die<num_rational::BigRational> {
    dice::from_fn(|mut fate| {
        let numerator = fate.roll(dice::u32(1..)).into();
        let denominator = fate.roll(dice::i32(1..)).into();
        num_rational::BigRational::new(numerator, denominator)
    })
}

/// Generates an [`Infinitesimal`] with a single one in an arbitrary dimension and zeros
/// for all other dimensions.
///
/// If the dimension count is zero, an [`Infinitesimal`] without any elements will be generated.
pub fn infinitesimal_variable<N: Number, I: Infinitesimal<N>>(dim: usize) -> impl Die<I> {
    dice::from_fn(move |mut fate| {
        if dim == 0 {
            I::zeros(0)
        } else {
            let idx = fate.roll(dice::uni_usize(0..dim));
            I::one(idx, dim)
        }
    })
}

/// Generates a [`Infinitesimal`] with arbitrary elements using [`Infinitesimal::from_dense`].
pub fn infinitesimal_dense<N: Number, I: Infinitesimal<N>, NDI: Die<N>>(
    dim: usize,
    infinitesimal_number_die: NDI,
) -> impl Die<I> {
    dice::from_fn(move |mut fate| {
        let elems = fate.roll(dice::vec(&infinitesimal_number_die, dim));
        I::from_dense(elems)
    })
}

/// Generates a [`Infinitesimal`] with arbitrary elements using [`Infinitesimal::from_sparse`].
pub fn infinitesimal_sparse<N: Number, I: Infinitesimal<N>, NDI: Die<N>>(
    dim: usize,
    infinitesimal_number_die: NDI,
) -> impl Die<I> {
    dice::from_fn(move |mut fate| {
        let number_with_index_die =
            dice::zip().two(dice::uni_usize(0..dim), &infinitesimal_number_die);
        let elems = fate.roll(dice::vec(number_with_index_die, 0..=dim));
        I::from_sparse(elems, dim)
    })
}

/// Generates an [`Infinitesimal`] with arbitrary elements.
pub fn infinitesimal<N: Number, I: Infinitesimal<N>, NDI: Die<N>>(
    dim: usize,
    infinitesimal_number_die: NDI,
) -> impl Die<I> {
    dice::from_fn(move |mut fate| {
        let infinitesimal_die = dice::one_of_die_once().three(
            infinitesimal_variable(dim),
            infinitesimal_dense(dim, &infinitesimal_number_die),
            infinitesimal_sparse(dim, &infinitesimal_number_die),
        );
        fate.roll(infinitesimal_die)
    })
}

/// Generates a [`Jet`] using [`Jet::constant`].
pub fn jet_constant<N: Number, I: Infinitesimal<N>, NDR: Die<N>>(
    dim: usize,
    real_number_die: NDR,
) -> impl Die<Jet<N, I>> {
    dice::from_fn(move |mut fate| {
        let real = fate.roll(&real_number_die);
        Jet::constant(real, dim)
    })
}

/// Generates a [`Jet`] using [`Jet::variable`].
pub fn jet_variable<N: Number, I: Infinitesimal<N>, NDR: Die<N>>(
    dim: usize,
    real_number_die: NDR,
) -> impl Die<Jet<N, I>> {
    dice::from_fn(move |mut fate| {
        let real = fate.roll(&real_number_die);
        let idx = fate.roll(dice::uni_usize(0..dim));
        Jet::variable(real, idx, dim)
    })
}

/// Generates a [`Jet`] with arbitrary values for the real part and the infinitesimal part.
pub fn jet<N: Number, I: Infinitesimal<N>, NDR: Die<N>, NDI: Die<N>>(
    dim: usize,
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

        let jet_die = dice::one_of_die_once().three(
            jet_constant(dim, &real_number_die),
            jet_variable(dim, &real_number_die),
            any_jet_die,
        );
        fate.roll(jet_die)
    })
}
