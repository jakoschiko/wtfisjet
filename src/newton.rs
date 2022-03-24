/// Provides an implementation for [Newton's method] that uses jets for calculating the derivative.
///
/// [Newton's method]: https://en.wikipedia.org/wiki/Newton%27s_method
use std::any::Any;

use num_traits::Float;
use rulinalg::{matrix::Matrix, vector::Vector};

use crate::{Dim, Infinitesimal, Jet, Number};

/// The multidimensional function whose root can be approximated with [`NewtonsMethod::find_root`].
///
/// You can interpret this type as a vector-function that has the same dimensionality for its
/// input vector and result vector.
///
/// The first argument is the dimensionality. It represents the actual length of the input and
/// the expected length of the result. Additionally you need to use this value if you create
/// jets inside the function.
///
/// The second argument is the input of the function. It represents the current guess of the
/// root finder. Use it to calculate the result.
///
/// The third argument is a buffer for the result of the function. The buffer is empty when
/// your function is called. Your function must push enough values into the buffer.
pub type MultidimensionalFunction<'a, N, I> =
    &'a mut dyn FnMut(Dim, &[Jet<N, I>], &mut Vec<Jet<N, I>>);

/// The result of [`NewtonsMethod::find_root`].
#[derive(Debug, Clone)]
pub struct NewtonsMethodResult<N: Number + Float + Any, I: Infinitesimal<N>> {
    /// Reason why the algorithm has stopped.
    pub state: NewtonsMethodState,
    /// The number of steps that were performed by the algorithm.
    pub steps: u32,
    /// The latest guess that were passed to the function.
    pub latest_guess: Vec<N>,
    /// The latest result that were produced by the function.
    pub latest_result: Vec<Jet<N, I>>,
}

/// Explains why the algorithm has stopped.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum NewtonsMethodState {
    /// The latest result has the wrong length.
    ///
    /// The length of the result produced by the function must match the provided [`Dim`].
    ResultHasWrongLength { wrong_len: usize },
    /// The latest result contains a jet with wrong dimension.
    ///
    /// Use only the provided [`Dim`] for constructing new jets inside the function.
    ResultHasJetWithWrongDim { wrong_dim: Dim },
    /// The latest result of the function produces an LES that was not solvable.
    ///
    /// This probably occurs because the function has no root.
    ResultNotSolvable,
    /// The maximum number of steps was performed without reaching the requested threshold.
    ///
    /// This probably occurs because the root finder found a minima. Please try to use another
    /// initial guess or check if your function is really continuously differentiable. Sometimes
    /// this is caused by floating-point errors. Then you could try to reduce the threshold
    /// or use a floating-point type with higher precision.
    MaxStepReached,
    /// The desired threshold was reached.
    ///
    /// This means that the root-finder has successfully found a root.
    ThresholdReached,
}

/// A multivariate root-finder based on [Newton's method].
///
/// It's multivariate because it operates on multiple variables at once.
///
/// [Newton's method]: https://en.wikipedia.org/wiki/Newton%27s_method
#[derive(Debug, Clone)]
pub struct NewtonsMethod<N: Number + Float + Any> {
    /// The maximum number of steps that will be performed.
    pub max_steps: u32,
    /// The threshold for deciding whether the result is close enough to zero.
    ///
    /// The algorithm applies the [maximum norm] to the result and then compares the output with
    /// the threshold. If the output is equal to or smaller than the threshold, the algorithm
    /// will stop.
    ///
    /// [maximum norm]: https://en.wikipedia.org/wiki/Norm_(mathematics)#Maximum_norm_(special_case_of:_infinity_norm,_uniform_norm,_or_supremum_norm)
    pub threshold: N,
}

impl<N: Number + Float + Any> NewtonsMethod<N> {
    /// Tries to approximate a root for the given functions.
    ///
    /// The algorithm will start with the given initial guess. Its length determines the
    /// dimensionality of this algorithm. The input passed to the function, the result of the
    /// function and the final result of the root-finding will all have the same length.
    ///
    /// The following guesses are calculated based on the result of the function and
    /// its derivative. The derivative is available because the functions takes and returns
    /// jets.
    pub fn find_root<I: Infinitesimal<N>>(
        &self,
        initial_guess: Vec<N>,
        function: MultidimensionalFunction<N, I>,
    ) -> NewtonsMethodResult<N, I> {
        let dim = Dim(initial_guess.len());
        let mut current_guess = initial_guess;
        let mut result = Vec::with_capacity(dim.0);
        let mut input = Vec::with_capacity(dim.0);
        let mut steps = 0;

        loop {
            // Initialize input
            input.clear();
            input.extend(
                current_guess
                    .iter()
                    .enumerate()
                    .map(|(idx, n)| Jet::variable(*n, idx, dim)),
            );

            // Evaluate the function
            result.clear();
            function(dim, &input, &mut result);

            // Check if result is valid
            if dim.0 != result.len() {
                return NewtonsMethodResult {
                    state: NewtonsMethodState::ResultHasWrongLength {
                        wrong_len: result.len(),
                    },
                    steps,
                    latest_guess: current_guess,
                    latest_result: result,
                };
            }
            for e in result.iter() {
                if dim != e.infinitesimal.dim() {
                    return NewtonsMethodResult {
                        state: NewtonsMethodState::ResultHasJetWithWrongDim {
                            wrong_dim: e.infinitesimal.dim(),
                        },
                        steps,
                        latest_guess: current_guess,
                        latest_result: result,
                    };
                }
            }

            // Check if result is close enough to zero
            let max_value = max_norm(result.iter().map(|e| &e.real));
            let threshold_reached = max_value <= self.threshold;

            if threshold_reached {
                return NewtonsMethodResult {
                    state: NewtonsMethodState::ThresholdReached,
                    steps,
                    latest_guess: current_guess,
                    latest_result: result,
                };
            }

            // Check if the maximal number of steps is reached
            if steps >= self.max_steps {
                return NewtonsMethodResult {
                    state: NewtonsMethodState::MaxStepReached,
                    steps,
                    latest_guess: current_guess,
                    latest_result: result,
                };
            }

            // Calculate the next guess
            if !next_guess(dim, &mut current_guess, &result) {
                return NewtonsMethodResult {
                    state: NewtonsMethodState::ResultNotSolvable,
                    steps,
                    latest_guess: current_guess,
                    latest_result: result,
                };
            }
            steps += 1;
        }
    }
}

fn max_norm<'a, N: 'a + Number, E: Iterator<Item = &'a N>>(elems: E) -> N {
    elems
        .map(N::abs)
        .max_by(N::total_cmp)
        .unwrap_or_else(N::zero)
}

fn next_guess<N: Number + Float + Any, I: Infinitesimal<N>>(
    dim: Dim,
    guess: &mut Vec<N>,
    result: &[Jet<N, I>],
) -> bool {
    let mut neg_f_x0_data = Vec::with_capacity(dim.0);
    let mut jacobi_data = Vec::with_capacity(dim.0 * dim.0);

    for result in result {
        neg_f_x0_data.push(-result.real);
        jacobi_data.extend(result.infinitesimal.dense_elems().cloned());
    }

    let neg_f_x0 = Vector::<N>::new(neg_f_x0_data);
    let jacobi = Matrix::<N>::new(dim.0, dim.0, jacobi_data);

    match jacobi.solve(neg_f_x0) {
        Err(_error) => false,
        Ok(x1_sub_x0) => {
            let x0 = Vector::<N>::from(guess.clone());
            let x1 = x1_sub_x0 + x0;

            guess.clear();
            guess.extend(x1.iter());
            true
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        newton::{NewtonsMethod, NewtonsMethodState},
        ConstInfinitesimal,
    };

    #[test]
    fn max_norm_examples() {
        assert_eq!(super::max_norm::<f32, _>([].iter()), 0.0);
        assert_eq!(super::max_norm::<f32, _>([1.0].iter()), 1.0);
        assert_eq!(super::max_norm::<f32, _>([-1.0].iter()), 1.0);
        assert_eq!(super::max_norm::<f32, _>([2.0, -1.0].iter()), 2.0);
        assert_eq!(super::max_norm::<f32, _>([-1.0, 2.0].iter()), 2.0);
        assert!(super::max_norm::<f32, _>([-1.0, f32::NAN].iter()).is_nan());
        assert!(super::max_norm::<f32, _>([f32::NAN, -1.0].iter()).is_nan());
    }

    // Approximates the supergolden ratio.
    //
    // Reference: https://personal.math.ubc.ca/~pwalls/math-python/roots-optimization/newton/#supergolden-ratio
    #[test]
    fn supergolden_ratio() {
        let newtons_method = NewtonsMethod {
            max_steps: 10,
            threshold: 1e-10,
        };

        let result = newtons_method.find_root::<ConstInfinitesimal<f64, 1>>(
            vec![1.0],
            &mut |_dim, input, result| {
                let x = &input[0];
                let r = x.clone().cube() - x.clone().square() - 1.0;
                result.push(r);
            },
        );

        println!("Result: {:#?}", result);

        assert_eq!(result.state, NewtonsMethodState::ThresholdReached);
        assert_eq!(result.steps, 6);
    }

    // Proves that the root-finder diverges if it tries to find a root for the root function.
    //
    // Reference: https://personal.math.ubc.ca/~pwalls/math-python/roots-optimization/newton/#divergent-example
    #[test]
    fn root_of_root_function() {
        let newtons_method = NewtonsMethod {
            max_steps: 100,
            threshold: 1e-2,
        };

        let result = newtons_method.find_root::<ConstInfinitesimal<f64, 1>>(
            vec![0.1],
            &mut |_dim, input, result| {
                let x = &input[0];
                let r = x.clone().pow_real(1.0 / 3.0);
                result.push(r);
            },
        );

        println!("Result: {:#?}", result);

        assert_eq!(result.state, NewtonsMethodState::MaxStepReached);
        assert_eq!(result.steps, newtons_method.max_steps);
    }
}
