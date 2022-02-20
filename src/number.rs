use std::fmt::Debug;
use std::ops::Neg;

/// The number type used by [`Jet`](crate::Jet).
///
/// Usually it's float like type, e.g. `f32`.
pub trait Number:
    num_traits::NumAssign
    + num_traits::Inv<Output = Self>
    + Neg<Output = Self>
    + Debug
    + Clone
    + PartialOrd
{
}

impl Number for f32 {}

impl Number for f64 {}

#[cfg(any(test, feature = "big-rational"))]
impl Number for num_rational::BigRational {}
