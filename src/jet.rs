#![allow(clippy::should_implement_trait)]

use num_traits::Inv;
use std::{
    fmt::Display,
    ops::{Add, Div, Mul, Neg, Sub},
};

// TODO: Impl assign operators

use crate::{Dim, Infinitesimal, Number};

/// Implements the jet, a n-dimensional dual number.
///
/// The jet consists of:
/// - A real part that represents the result of a function.
/// - An infinitesimal part that represents the exact derivative for each input of the function.
///
/// # Using a jet based function
///
/// If you want to use a jet based function and calculate the derivatives with respect to
/// some of its inputs, you probably want to use [`Jet::variable`] to initialize these inputs.
/// Use a different index for each input.
///
/// # Writing a jet based function
///
/// If you want to write a jet based function, you need to be careful. There are some pitfalls:
/// - Jets do **not** fulfill the properties of a mathematical field. Most notably there are
/// multiple jets without an inverse (all jets with zero in their real part). You need to
/// consider this if you divide jets.
/// - Continuously differentiable functions are very nice for algorithms like Newton's method.
/// Unfortunately you can easily lose this property, e.g. by using if/else.
/// - There is no obvious order for jets.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Jet<N: Number, I: Infinitesimal<N>> {
    /// The real part of the jet.
    pub real: N,
    /// The infinitesimal part of the jet.
    pub infinitesimal: I,
}

impl<N: Number, I: Infinitesimal<N>> Jet<N, I> {
    /// Creates an instance with the given real and infinitesimal parts.
    pub fn new(real: N, infinitesimal: I) -> Self {
        Self {
            real,
            infinitesimal,
        }
    }

    /// Creates an instance with the given real part and zeros for the infinitesimal part.
    ///
    /// This represents a jet that has no dependencies to any variables.
    ///
    /// # Panic
    ///
    /// This function panics if the implementation for `Infinitesimal` does no support the
    /// dimension count.
    pub fn constant(real: N, dim: Dim) -> Self {
        Self {
            real,
            infinitesimal: I::zeros(dim),
        }
    }

    /// Creates an instance with the given real part and a single one for the given dimension
    /// in the infinitesimal part.
    ///
    /// Useful for initializing the inputs of the function for which you want to calculate the
    /// derivatives. Use a different index for each input.
    ///
    /// # Panic
    ///
    /// This function panics if
    /// - the given index is equal to or greater than the dimension count or
    /// - the implementation for `Infinitesimal` does not support the dimension count.
    pub fn variable(real: N, idx: usize, dim: Dim) -> Self {
        Self {
            real,
            infinitesimal: I::one(idx, dim),
        }
    }

    /// Adds both jets together.
    pub fn add(self, rhs: Self) -> Self {
        Self {
            real: self.real + rhs.real,
            infinitesimal: self.infinitesimal.add(rhs.infinitesimal),
        }
    }

    /// Adds the real number to the jet.
    pub fn add_real(self, rhs: N) -> Self {
        Self {
            real: self.real + rhs,
            infinitesimal: self.infinitesimal,
        }
    }

    /// Subtracts the right jet from the left jet.
    pub fn sub(self, rhs: Self) -> Self {
        Self {
            real: self.real - rhs.real,
            infinitesimal: self.infinitesimal.sub(rhs.infinitesimal),
        }
    }

    /// Subtracts the real number from the jet.
    pub fn sub_real(self, rhs: N) -> Self {
        Self {
            real: self.real - rhs,
            infinitesimal: self.infinitesimal,
        }
    }

    /// Negates the jet.
    pub fn neg(self) -> Self {
        Self {
            real: -self.real,
            infinitesimal: self.infinitesimal.neg(),
        }
    }

    /// Multiplies both jets.
    pub fn mul(self, rhs: Self) -> Self {
        let infinitesimal = self
            .infinitesimal
            .scale(rhs.real.clone())
            .add(rhs.infinitesimal.scale(self.real.clone()));

        Self {
            real: self.real * rhs.real,
            infinitesimal,
        }
    }

    /// Multiplies the real number with the jet.
    pub fn mul_real(self, rhs: N) -> Self {
        let infinitesimal = self.infinitesimal.scale(rhs.clone());

        Self {
            real: self.real * rhs,
            infinitesimal,
        }
    }

    /// Divides the left jet by the right jet if the right jet is invertible.
    ///
    /// The right jet is not invertible if its **real part** is zero even if its infinitesimal
    /// part is non-zero. I.e. there are multiple jets without an inverse. Therefore jets do
    /// **not** fulfill the properties of a mathematical field.
    pub fn div(self, rhs: Self) -> Self {
        let lhs_inv = rhs.real.inv();
        let real = self.real * lhs_inv.clone();
        let infinitesimal = self
            .infinitesimal
            .sub(rhs.infinitesimal.scale(real.clone()))
            .scale(lhs_inv);

        Self {
            real,
            infinitesimal,
        }
    }

    /// Divides the jet by the real number if the real number is invertible.
    ///
    /// The real number is not invertible if it is zero.
    pub fn div_real(self, rhs: N) -> Self {
        let lhs_inv = rhs.inv();
        let real = self.real * lhs_inv.clone();
        let infinitesimal = self.infinitesimal.scale(lhs_inv);

        Self {
            real,
            infinitesimal,
        }
    }

    /// Divides the left jet by the right jet if the right jet is invertible or returns `None`
    /// otherwise.
    ///
    /// The right jet is not invertible if its **real part** is zero even if its infinitesimal
    /// part is non-zero. I.e. there are multiple jets without an inverse. Therefore jets do
    /// **not** fulfill the properties of a mathematical field.
    pub fn checked_div(self, rhs: Self) -> Option<Self> {
        if rhs.real.is_zero() {
            None
        } else {
            Some(self.div(rhs))
        }
    }

    /// Divides the jet by the real number if the real number is invertible or returns `None`
    /// otherwise.
    ///
    /// The real number is not invertible if it is zero.
    pub fn checked_div_real(self, rhs: N) -> Option<Self> {
        if rhs.is_zero() {
            None
        } else {
            Some(self.div_real(rhs))
        }
    }

    /// Returns the inverse of the jet if it is invertible.
    ///
    /// The jet is not invertible if its **real part** is zero even if its infinitesimal
    /// part is non-zero. I.e. there are multiple jets without an inverse. Therefore jets do
    /// **not** fulfill the properties of a mathematical field.
    pub fn inv(self) -> Self {
        let inv = self.real.inv();
        let infinitesimal = self.infinitesimal.scale(-inv.clone() * inv.clone());

        Self {
            real: inv,
            infinitesimal,
        }
    }

    /// Returns the inverse of the jet if it is invertible or `None` otherwise.
    ///
    /// The jet is not invertible if its **real part** is zero even if its infinitesimal
    /// part is non-zero. I.e. there are multiple jets without an inverse. Therefore jets do
    /// **not** fulfill the properties of a mathematical field.
    pub fn checked_inv(self) -> Option<Self> {
        if self.real.is_zero() {
            None
        } else {
            Some(self.inv())
        }
    }
}

impl<N: Number, I: Infinitesimal<N>> Add<Self> for Jet<N, I> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        self.add(rhs)
    }
}

impl<N: Number, I: Infinitesimal<N>> Add<N> for Jet<N, I> {
    type Output = Self;

    fn add(self, rhs: N) -> Self::Output {
        self.add_real(rhs)
    }
}

impl<N: Number, I: Infinitesimal<N>> Sub<Self> for Jet<N, I> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        self.sub(rhs)
    }
}

impl<N: Number, I: Infinitesimal<N>> Sub<N> for Jet<N, I> {
    type Output = Self;

    fn sub(self, rhs: N) -> Self::Output {
        self.sub_real(rhs)
    }
}

impl<N: Number, I: Infinitesimal<N>> Neg for Jet<N, I> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        self.neg()
    }
}

impl<N: Number, I: Infinitesimal<N>> Mul<Self> for Jet<N, I> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        self.mul(rhs)
    }
}

impl<N: Number, I: Infinitesimal<N>> Mul<N> for Jet<N, I> {
    type Output = Self;

    fn mul(self, rhs: N) -> Self::Output {
        self.mul_real(rhs)
    }
}

impl<N: Number, I: Infinitesimal<N>> Div<Self> for Jet<N, I> {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        self.div(rhs)
    }
}

impl<N: Number, I: Infinitesimal<N>> Div<N> for Jet<N, I> {
    type Output = Self;

    fn div(self, rhs: N) -> Self::Output {
        self.div_real(rhs)
    }
}

impl<N: Number, I: Infinitesimal<N>> Inv for Jet<N, I> {
    type Output = Self;

    fn inv(self) -> Self::Output {
        self.inv()
    }
}

impl<N: Number + Display, I: Infinitesimal<N>> Display for Jet<N, I> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({} + [", self.real)?;

        let dim = self.infinitesimal.dim();
        for (idx, elem) in self.infinitesimal.dense_elems().enumerate() {
            if idx + 1 < dim.0 {
                write!(f, "{}, ", elem)?;
            } else {
                write!(f, "{}", elem)?;
            }
        }

        write!(f, "]h)")
    }
}

#[cfg(test)]
mod tests {
    use diceprop::{props, Elem, Fun1, Fun2, Set};
    use dicetest::prelude::*;
    use num_rational::BigRational;
    use num_traits::{Inv, One, Zero};

    use crate::{DenseInfinitesimal, Dim, Infinitesimal, Jet};

    type TestInfinitesimal = DenseInfinitesimal<BigRational>;

    type TestJet = Jet<BigRational, TestInfinitesimal>;

    fn rational(num: i32, denom: u32) -> BigRational {
        BigRational::new(num.into(), denom.into())
    }

    fn infinitesimal<const D: usize>(elems: [BigRational; D]) -> TestInfinitesimal {
        TestInfinitesimal::from_dense(elems)
    }

    fn jet<const D: usize>(
        real_part: BigRational,
        infinitesimal_part: [BigRational; D],
    ) -> TestJet {
        TestJet::new(real_part, infinitesimal(infinitesimal_part))
    }

    fn number_die() -> impl Die<BigRational> {
        crate::dice::big_rational_number()
    }

    fn non_zero_number_die() -> impl Die<BigRational> {
        crate::dice::big_rational_non_zero_number()
    }

    fn jet_die(dim: Dim) -> impl Die<TestJet> {
        crate::dice::jet(dim, number_die(), number_die())
    }

    fn non_zero_jet_die(dim: Dim) -> impl Die<TestJet> {
        crate::dice::jet(dim, non_zero_number_die(), number_die())
    }

    #[test]
    fn add_examples() {
        let x = jet(rational(2, 1), [rational(3, 1), rational(0, 1)]);
        let y = jet(rational(4, 1), [rational(0, 1), rational(2, 1)]);
        let z = jet(rational(6, 1), [rational(3, 1), rational(2, 1)]);
        assert_eq!(x + y, z);
    }

    #[test]
    fn add_real_examples() {
        let x = jet(rational(2, 1), [rational(3, 1), rational(0, 1)]);
        let y = rational(4, 1);
        let z = jet(rational(6, 1), [rational(3, 1), rational(0, 1)]);
        assert_eq!(x + y, z);
    }

    #[test]
    fn sub_examples() {
        let x = jet(rational(2, 1), [rational(3, 1), rational(0, 1)]);
        let y = jet(rational(4, 1), [rational(0, 1), rational(2, 1)]);
        let z = jet(rational(-2, 1), [rational(3, 1), rational(-2, 1)]);
        assert_eq!(x - y, z);
    }

    #[test]
    fn sub_real_examples() {
        let x = jet(rational(2, 1), [rational(3, 1), rational(0, 1)]);
        let y = rational(4, 1);
        let z = jet(rational(-2, 1), [rational(3, 1), rational(0, 1)]);
        assert_eq!(x - y, z);
    }

    #[test]
    fn neg_examples() {
        let x = jet(rational(2, 1), [rational(3, 1), rational(0, 1)]);
        let y = jet(rational(-2, 1), [rational(-3, 1), rational(0, 1)]);
        assert_eq!(-x, y);
    }

    #[test]
    fn mul_examples() {
        let x = jet(rational(2, 1), [rational(3, 1), rational(0, 1)]);
        let y = jet(rational(4, 1), [rational(0, 1), rational(2, 1)]);
        let z = jet(rational(8, 1), [rational(12, 1), rational(4, 1)]);
        assert_eq!(x * y, z);
    }

    #[test]
    fn mul_real_examples() {
        let x = jet(rational(2, 1), [rational(3, 1), rational(0, 1)]);
        let y = rational(4, 1);
        let z = jet(rational(8, 1), [rational(12, 1), rational(0, 1)]);
        assert_eq!(x * y, z);
    }

    #[test]
    fn div_examples() {
        let x = jet(rational(2, 1), [rational(3, 1), rational(0, 1)]);
        let y = jet(rational(4, 1), [rational(0, 1), rational(2, 1)]);
        let z = jet(rational(1, 2), [rational(3, 4), rational(-1, 4)]);
        assert_eq!(x / y, z);
    }

    #[test]
    fn div_real_examples() {
        let x = jet(rational(2, 1), [rational(3, 1), rational(0, 1)]);
        let y = rational(4, 1);
        let z = jet(rational(1, 2), [rational(3, 4), rational(0, 1)]);
        assert_eq!(x / y, z);
    }

    #[test]
    fn checked_div_examples() {
        {
            let x = jet(rational(2, 1), [rational(3, 1), rational(0, 1)]);
            let y = jet(rational(4, 1), [rational(0, 1), rational(2, 1)]);
            let z = jet(rational(1, 2), [rational(3, 4), rational(-1, 4)]);
            assert_eq!(Jet::checked_div(x, y), Some(z));
        }
        {
            let x = jet(rational(2, 1), [rational(3, 1), rational(0, 1)]);
            let y = jet(rational(0, 1), [rational(0, 1), rational(2, 1)]);
            assert_eq!(Jet::checked_div(x, y), None);
        }
    }

    #[test]
    fn checked_div_real_examples() {
        {
            let x = jet(rational(2, 1), [rational(3, 1), rational(0, 1)]);
            let y = rational(4, 1);
            let z = jet(rational(1, 2), [rational(3, 4), rational(0, 1)]);
            assert_eq!(Jet::checked_div_real(x, y), Some(z));
        }
        {
            let x = jet(rational(2, 1), [rational(3, 1), rational(0, 1)]);
            let y = rational(0, 1);
            assert_eq!(Jet::checked_div_real(x, y), None);
        }
    }

    #[test]
    fn inv_examples() {
        let x = jet(rational(2, 1), [rational(3, 1), rational(0, 1)]);
        let y = jet(rational(1, 2), [rational(-3, 4), rational(0, 4)]);
        assert_eq!(<TestJet as Inv>::inv(x), y);
    }

    #[test]
    fn checked_inv_examples() {
        {
            let x = jet(rational(2, 1), [rational(3, 1), rational(0, 1)]);
            let y = jet(rational(1, 2), [rational(-3, 4), rational(0, 4)]);
            assert_eq!(Jet::checked_inv(x), Some(y));
        }
        {
            let x = jet(rational(0, 1), [rational(3, 1), rational(0, 1)]);
            assert_eq!(Jet::checked_inv(x), None);
        }
    }

    #[test]
    fn display_examples() {
        type MyJet = Jet<f32, DenseInfinitesimal<f32>>;

        assert_eq!(
            "(2.4 + []h)",
            format!("{}", MyJet::new(2.4, DenseInfinitesimal::from_dense([]))),
        );

        assert_eq!(
            "(-0.2 + [4]h)",
            format!("{}", MyJet::new(-0.2, DenseInfinitesimal::from_dense([4.]))),
        );

        assert_eq!(
            "(10 + [0, -0.1]h)",
            format!(
                "{}",
                MyJet::new(10., DenseInfinitesimal::from_dense([0., -0.1]))
            ),
        );

        assert_eq!(
            "(0 + [-1, 0.5, 0]h)",
            format!(
                "{}",
                MyJet::new(0., DenseInfinitesimal::from_dense([-1., 0.5, 0.]))
            ),
        );
    }

    #[test]
    fn commutative_ring() {
        Dicetest::repeatedly().run(|mut fate| {
            let dim = fate.roll(crate::dice::dim());

            let set = Set::new("Jet", jet_die(dim));
            let vars = fate.roll(set.vars(["x", "y", "z"]));
            let add = Fun2::new("Jet::add", Jet::add);
            let mul = Fun2::new("Jet::mul", Jet::mul);
            let neg = Fun1::new("Jet::neg", Jet::neg);
            let zero = Elem::new("zero", Jet::constant(BigRational::zero(), dim));
            let one = Elem::new("one", Jet::constant(BigRational::one(), dim));

            props::algebra::commutative_ring(vars, add, mul, neg, zero, one);
        });
    }

    #[test]
    fn add_real_is_equal_to_add_with_constant() {
        Dicetest::repeatedly().run(|mut fate| {
            let dim = fate.roll(crate::dice::dim());

            let x = fate.roll(jet_die(dim));
            let a = fate.roll(number_die());

            assert_eq!(
                Jet::add_real(x.clone(), a.clone()),
                Jet::add(x, Jet::constant(a, dim))
            );
        });
    }

    #[test]
    fn sub_is_equal_to_add_with_neg() {
        Dicetest::repeatedly().run(|mut fate| {
            let dim = fate.roll(crate::dice::dim());

            let [x, y] = fate.roll(dice::array(jet_die(dim)));

            assert_eq!(Jet::sub(x.clone(), y.clone()), Jet::add(x, Jet::neg(y)));
        });
    }

    #[test]
    fn sub_real_is_equal_to_sub_with_constant() {
        Dicetest::repeatedly().run(|mut fate| {
            let dim = fate.roll(crate::dice::dim());

            let x = fate.roll(jet_die(dim));
            let a = fate.roll(number_die());

            assert_eq!(
                Jet::sub_real(x.clone(), a.clone()),
                Jet::sub(x, Jet::constant(a, dim))
            );
        });
    }

    #[test]
    fn mul_real_is_equal_to_mul_with_constant() {
        Dicetest::repeatedly().run(|mut fate| {
            let dim = fate.roll(crate::dice::dim());

            let x = fate.roll(jet_die(dim));
            let a = fate.roll(number_die());

            assert_eq!(
                Jet::mul_real(x.clone(), a.clone()),
                Jet::mul(x, Jet::constant(a, dim))
            );
        });
    }

    #[test]
    fn div_is_right_inverse_of_mul() {
        Dicetest::repeatedly().run(|mut fate| {
            let dim = fate.roll(crate::dice::dim());

            let x = fate.roll(jet_die(dim));
            let y = fate.roll(non_zero_jet_die(dim));

            let xy = Jet::mul(x.clone(), y.clone());

            assert_eq!(x, Jet::div(xy, y));
        });
    }

    #[test]
    fn div_real_is_equal_to_div_with_constant() {
        Dicetest::repeatedly().run(|mut fate| {
            let dim = fate.roll(crate::dice::dim());

            let x = fate.roll(jet_die(dim));
            let a = fate.roll(non_zero_number_die());

            assert_eq!(
                Jet::div_real(x.clone(), a.clone()),
                Jet::div(x, Jet::constant(a, dim))
            );
        });
    }

    #[test]
    fn checked_div_real_is_equal_to_checked_div_with_constant() {
        Dicetest::repeatedly().run(|mut fate| {
            let dim = fate.roll(crate::dice::dim());

            let x = fate.roll(jet_die(dim));
            let a = fate.roll(number_die());

            assert_eq!(
                Jet::checked_div_real(x.clone(), a.clone()),
                Jet::checked_div(x, Jet::constant(a, dim))
            );
        });
    }

    #[test]
    fn div_is_equal_to_mul_with_inv() {
        Dicetest::repeatedly().run(|mut fate| {
            let dim = fate.roll(crate::dice::dim());

            let x = fate.roll(jet_die(dim));
            let y = fate.roll(non_zero_jet_die(dim));

            assert_eq!(Jet::div(x.clone(), y.clone()), Jet::mul(x, Jet::inv(y)));
        });
    }
}
