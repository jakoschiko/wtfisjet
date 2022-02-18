/// The number type used by [`Jet`](crate::Jet).
///
/// Usually it's float like type, e.g. `f32`.
pub trait Number: num_traits::NumAssign + num_traits::Inv<Output = Self> + PartialOrd {}

impl Number for f32 {}

impl Number for f64 {}

#[cfg(feature = "num-rational")]
impl Number for num_rational::BigRational {}
