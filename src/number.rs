use std::cmp::Ordering;
use std::fmt::Debug;
use std::ops::Neg;

/// The number type used by [`Jet`](crate::Jet).
///
/// This library assumes that the number type has at least [`field like properties`].
/// Some parts of library have even stricter requirements.
///
/// For production use you probably want to use either `f32` or `f64`, even if they are
/// technically not a field. These types are also the reason why this trait uses the trait bound
/// `PartialOrd` and not `Ord`. For testing purposes you might want to use
/// `num_rational::BigRational` because it has arbitrary precision which relieves you from the
/// burden of dealing with epsilons.
///
/// [`field like properties`]: https://en.wikipedia.org/wiki/Field_(mathematics)#Classic_definition
pub trait Number:
    num_traits::NumAssign
    + num_traits::Inv<Output = Self>
    + num_traits::Signed
    + Neg<Output = Self>
    + Debug
    + Clone
    + PartialOrd
{
    fn two() -> Self;

    fn three() -> Self;

    fn from_integer(integer: i32) -> Self;

    /// Compares both numbers using a total order.
    ///
    /// This function must be consistent with the number's [`PartialOrd`] implementation.
    /// Unfortunately we can't require [`Ord`] because it's not implemented for the
    /// floating point number types [`f32`] and [`f64`]. The usage of [`PartialOrd`] should
    /// be preferred except you really need a total order.
    fn total_cmp(&self, rhs: &Self) -> Ordering;

    /// Same as [`Neg::neg`](std::ops::Neg::neg), but in place.
    fn neg_in_place(&mut self);
}

impl Number for f32 {
    #[inline]
    fn two() -> Self {
        2.0
    }

    #[inline]
    fn three() -> Self {
        3.0
    }

    #[inline]
    fn from_integer(integer: i32) -> Self {
        integer as Self
    }

    #[inline]
    fn total_cmp(&self, rhs: &Self) -> Ordering {
        f32::total_cmp(self, rhs)
    }

    #[inline]
    fn neg_in_place(&mut self) {
        *self = -*self;
    }
}

impl Number for f64 {
    #[inline]
    fn two() -> Self {
        2.0
    }

    #[inline]
    fn three() -> Self {
        3.0
    }

    #[inline]
    fn from_integer(integer: i32) -> Self {
        integer as Self
    }

    #[inline]
    fn total_cmp(&self, rhs: &Self) -> Ordering {
        f64::total_cmp(self, rhs)
    }

    #[inline]
    fn neg_in_place(&mut self) {
        *self = -*self;
    }
}

#[cfg(any(test, feature = "big-rational-number"))]
impl Number for num_rational::BigRational {
    fn two() -> Self {
        Self::from_integer(2.into())
    }

    fn three() -> Self {
        Self::from_integer(3.into())
    }

    fn from_integer(integer: i32) -> Self {
        Self::from_integer(integer.into())
    }

    fn total_cmp(&self, rhs: &Self) -> Ordering {
        Ord::cmp(self, rhs)
    }

    fn neg_in_place(&mut self) {
        use num_traits::Zero;

        let temp = std::mem::replace(self, Self::zero());
        *self = -temp;
    }
}
