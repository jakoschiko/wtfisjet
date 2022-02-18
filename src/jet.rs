use crate::{Infinitesimal, Number};

/// Implements the jet, a n-dimensional dual number.
///
/// The jet consists of:
/// - A real part that represents the result of a function.
/// - An infinitesimal part that represents the exact derivative for each input of the function.
///
/// If you want to write functions with jets, you need to be careful. There are some pitfalls:
/// - Jets do **not** fulfil the properties of a mathematical field. Most notably there are
/// multiple jets without an inverse (all jets with zero in their real part). You need to
/// consider this if you divide jets.
/// - There is no obvious order for jets.
/// - The derivatives are only plausible if your function is continuously differentiable.
/// You can easily lose this property, e.g. by using if/else.
#[derive(Debug, Clone, PartialEq)]
pub struct Jet<N: Number, I: Infinitesimal<N>> {
    /// The real part of the jet.
    pub real: N,
    /// The infinitesimal part of the jet.
    pub infinitesimal: I,
}

impl<N: Number, I: Infinitesimal<N>> Jet<N, I> {
    /// Creates a new instance with the given real and infinitesimal parts.
    pub fn new(real: N, infinitesimal: I) -> Self {
        Self {
            real,
            infinitesimal,
        }
    }
}
