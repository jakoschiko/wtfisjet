use crate::{Infinitesimal, Number};

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
/// - Jets do **not** fulfil the properties of a mathematical field. Most notably there are
/// multiple jets without an inverse (all jets with zero in their real part). You need to
/// consider this if you divide jets.
/// - Continuously differentiable functions are very nice for algorithms like Newton's method.
/// Unfortunately you can easily lose this property, e.g. by using if/else.
/// - There is no obvious order for jets.
#[derive(Debug, Clone, PartialEq)]
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
    pub fn constant(real: N, dim: usize) -> Self {
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
    pub fn variable(real: N, idx: usize, dim: usize) -> Self {
        Self {
            real,
            infinitesimal: I::one(idx, dim),
        }
    }
}
