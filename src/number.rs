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
    + Neg<Output = Self>
    + Debug
    + Clone
    + PartialOrd
{
    /// Same as [`Neg::neg`](std::ops::Neg::neg), but in place.
    fn neg_in_place(&mut self);
}

impl Number for f32 {
    #[inline]
    fn neg_in_place(&mut self) {
        *self = -*self;
    }
}

impl Number for f64 {
    #[inline]
    fn neg_in_place(&mut self) {
        *self = -*self;
    }
}

#[cfg(any(test, feature = "big-rational-number"))]
impl Number for num_rational::BigRational {
    fn neg_in_place(&mut self) {
        use num_traits::Zero;

        let temp = std::mem::replace(self, Self::zero());
        *self = -temp;
    }
}
