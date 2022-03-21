//! Provides an implementation for B-splines that uses jets as control points.
//!
//! Finding control points that fulfill the required properties is a non-trivial task.
//! One approach is to use a numeric solver to approximate the solution. With jets we
//! can calculate the derivatives with respect to the control points. This allows us to
//! use very efficient numeric solvers (e.g. Newton's method).

use std::{cmp::Ordering, fmt::Display, iter::repeat};

use crate::{Dim, Infinitesimal, Jet, Number};

/// Error that occurred during [`BSpline`] construction.
#[derive(Debug, Clone)]
pub struct BSplineError<N: Number> {
    /// The reason why the construction failed.
    pub reason: BSplineErrorReason,
    /// The requested degree of the [`BSpline`].
    pub degree: usize,
    /// The knots that were used to construct the [`BSpline`].
    pub knots: Vec<N>,
}

/// The reason why [`BSplineError`] was returned.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BSplineErrorReason {
    /// There were not enough knots for constructing the [`BSpline`] with the requested degree.
    NotEnoughKnots,
    /// The knots needs to be strictly increasing, but there is at least one knot that violates
    /// this property.
    NotStrictlyIncreasingKnots {
        /// The index of the knot that violates this property.
        wrong_knot_index: usize,
    },
}

impl<N: Number> Display for BSplineError<N> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.reason {
            BSplineErrorReason::NotEnoughKnots => {
                write!(
                    f,
                    "A B-spline with degree {} needs at least {} knots, but {} knots were given",
                    self.degree,
                    min_knot_count(self.degree),
                    self.knots.len(),
                )
            }
            BSplineErrorReason::NotStrictlyIncreasingKnots { wrong_knot_index } => {
                write!(
                    f,
                    "The knots of the B-spline must be strictly increasing, but the knot at index {} is either smaller than or equal to its predecessor",
                    wrong_knot_index,
                )
            }
        }
    }
}

impl<N: Number, I: Infinitesimal<N>> Display for BSplineCurveError<N, I> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let knot_count = self.bspline.knot_count();
        write!(
            f,
            "A B-spline curve with {} knots needs exactly {} control points, but {} control points were given",
            knot_count,
            necessary_control_point_count(self.bspline.degree, knot_count),
            self.control_points.len(),
        )
    }
}

/// A [B-spline] that constructed from a set of knots.
///
/// Creating a [`BSpline`] is the first step of creating a [`BSplineCurve`]. The knots define
/// the borders of the polynomial function pieces. The next step is to choose control points
/// to construct a [`BSplineCurve`] that has the desired properties (interpolated points,
/// specific derivative, etc.). There are different strategies for choosing the right
/// control points, though this library encourage you to use a numeric solver to find a solution.
///
/// [B-spline]: https://en.wikipedia.org/wiki/B-spline
#[derive(Debug, Clone)]
pub struct BSpline<N: Number> {
    // The degree of the polynomial function pieces.
    degree: usize,
    // This consists of
    // - The original knots passed to the constructor and
    // - `degree` padding knots at the front and
    // - `degree` paddings knots at the end.
    // The paddings knots are necessary for De Boor's algorithm. We create them by simply
    // repeating the first and last original knot.
    padded_knots: Vec<N>,
}

impl<N: Number> BSpline<N> {
    /// Tries to create a B-spline.
    ///
    /// The degree parameters defines the maximal degree of the the polynomial function pieces
    /// and the knots define their borders. There must be at least `degree + 1` knots and
    /// they must be strictly increasing.
    pub fn new(degree: usize, knots: Vec<N>) -> Result<Self, BSplineError<N>> {
        let knots = knots.into_iter().collect::<Vec<_>>();

        if knots.len() < min_knot_count(degree) {
            return Err(BSplineError {
                reason: BSplineErrorReason::NotEnoughKnots,
                degree,
                knots,
            });
        }

        let strictly_increasing_violation = knots.windows(2).enumerate().find(|(_, ks)| {
            !matches!(
                PartialOrd::partial_cmp(&ks[0], &ks[1]),
                Some(Ordering::Less)
            )
        });

        if let Some((i, _)) = strictly_increasing_violation {
            return Err(BSplineError {
                reason: BSplineErrorReason::NotStrictlyIncreasingKnots {
                    wrong_knot_index: i + 1,
                },
                degree,
                knots,
            });
        }

        let padded_knots = if knots.is_empty() {
            knots
        } else {
            let first = knots.first().unwrap();
            let last = knots.last().unwrap();
            let padding_init = repeat(first.clone()).take(degree);
            let padding_tail = repeat(last.clone()).take(degree);
            padding_init.chain(knots).chain(padding_tail).collect()
        };

        Ok(Self {
            degree,
            padded_knots,
        })
    }

    /// The degree of the polynomial function pieces of the B-spline.
    pub fn degree(&self) -> usize {
        self.degree
    }

    /// Returns the number of knots that were used to create the B-spline.
    pub fn knot_count(&self) -> usize {
        self.padded_knots.len() - 2 * self.degree
    }

    /// The knots that were used to create the B-spline.
    pub fn knots(&self) -> &[N] {
        &self.padded_knots[self.degree..self.padded_knots.len() - self.degree]
    }

    fn find_interval(&self, x: &N) -> usize {
        let knots = self.knots();

        // Calculate the knots between the intervals
        let interval_knots = if self.degree == 0 {
            // Because `degree` is at least 0, we know that we have at least 1 knot
            debug_assert!(!knots.is_empty());

            // `degree` 0 is asymmetrical and "right-heavy" and we are not allowed to
            // omit the last knot, therefore omit only the first knot
            &knots[1..]
        } else {
            // Because `degree` is at least 1, we know that we have at least 2 knots
            debug_assert!(knots.len() >= 2);

            // We omit the first and last knot
            &knots[1..knots.len() - 1]
        };

        // Find the interval that contains `x`
        let interval_knot_index = match interval_knots.binary_search_by(|k| N::total_cmp(k, x)) {
            Ok(index) => index,
            Err(index) => index,
        };

        interval_knot_index + self.degree
    }

    fn value<I: Infinitesimal<N>>(
        &self,
        control_points: &[Jet<N, I>],
        buffer: &mut BSplineCurveBuffer<N, I>,
        x: &N,
    ) -> Jet<N, I> {
        // The algorithm is inspired by https://en.wikipedia.org/wiki/De_Boor%27s_algorithm

        // We assume that the number of control points correct because it was checked before
        debug_assert_eq!(
            control_points.len(),
            necessary_control_point_count(self.degree, self.knot_count()),
        );

        let interval_index = self.find_interval(x);

        let relevant_control_points = control_points[interval_index - self.degree..=interval_index]
            .iter()
            .cloned();

        buffer.0.clear();
        buffer.0.extend(relevant_control_points);

        for step in 1..=self.degree {
            for i in (step..=self.degree).rev() {
                let left_knot_index = interval_index + i - self.degree;
                let left_knot = &self.padded_knots[left_knot_index];

                let right_knot_index = interval_index + 1 + i - step;
                let right_knot = &self.padded_knots[right_knot_index];

                let knot_diff = right_knot.clone() - left_knot.clone();

                // `knot_diff` is non-zero because we assume that the knots are strictly increasing
                // and at most one padding knot is used for the difference
                debug_assert!(!knot_diff.is_zero());

                let alpha = (x.clone() - left_knot.clone()) / knot_diff;
                let beta = N::one() - alpha.clone();

                buffer.0[i] = buffer.0[i - 1].clone() * beta + buffer.0[i].clone() * alpha;
            }
        }

        buffer.0.pop().unwrap()
    }

    fn derivative<I: Infinitesimal<N>>(
        &self,
        dim: Dim,
        control_points: &[Jet<N, I>],
        buffer: &mut BSplineCurveBuffer<N, I>,
        x: &N,
    ) -> Jet<N, I> {
        // The algorithm is inspired by https://en.wikipedia.org/wiki/De_Boor%27s_algorithm

        // We assume that the number of control points correct because it was checked before
        debug_assert_eq!(
            control_points.len(),
            necessary_control_point_count(self.degree, self.knot_count()),
        );

        if self.degree == 0 {
            return Jet::new(N::zero(), I::zeros(dim));
        }

        let interval_index = self.find_interval(x);

        let derivative_control_points = (0..self.degree).map(|j| {
            let left_control_point = &control_points[j + interval_index - self.degree];
            let right_control_point = &control_points[j + interval_index - self.degree + 1];
            let control_point_diff = right_control_point.clone() - left_control_point.clone();

            let left_knot = &self.padded_knots[j + interval_index - self.degree + 1];
            let right_knot = &self.padded_knots[j + interval_index + 1];
            let knot_diff = right_knot.clone() - left_knot.clone();

            // `knot_diff` is non-zero because we assume that the knots are strictly increasing
            // and at most one padding knot is used for the difference
            debug_assert!(!knot_diff.is_zero());

            let degree_number = N::from_integer(self.degree as i32);

            control_point_diff * degree_number / knot_diff
        });

        buffer.0.clear();
        buffer.0.extend(derivative_control_points);

        for step in 1..self.degree {
            for i in (step..self.degree).rev() {
                let left_knot = &self.padded_knots[interval_index + i - (self.degree - 1)];
                let right_knot = &self.padded_knots[interval_index + 1 + i - step];
                let knot_diff = right_knot.clone() - left_knot.clone();

                // `knot_diff` is non-zero because we assume that the knots are strictly increasing
                // and at most one padding knot is used for the difference
                debug_assert!(!knot_diff.is_zero());

                let alpha = (x.clone() - left_knot.clone()) / knot_diff;
                let beta = N::one() - alpha.clone();

                buffer.0[i] = buffer.0[i - 1].clone() * beta + buffer.0[i].clone() * alpha;
            }
        }

        buffer.0.pop().unwrap()
    }

    /// The number of control points that are necessary for construction a [`BSplineCurve`].
    ///
    /// The necessary number of control points is `knot_count + max(0, degree - 1)`.
    pub fn necessary_control_point_count(&self) -> usize {
        necessary_control_point_count(self.degree, self.knot_count())
    }
}

/// Buffer for temporary values during the [`BSplineCurve`] evaluation.
///
/// Exists only for performance reasons.
#[derive(Debug, Clone)]
pub struct BSplineCurveBuffer<N: Number, I: Infinitesimal<N>>(Vec<Jet<N, I>>);

impl<N: Number, I: Infinitesimal<N>> BSplineCurveBuffer<N, I> {
    /// Creates a buffer with an ideal preallocated size for B-spline curves with the given degree.
    pub fn new(degree: usize) -> Self {
        Self(Vec::with_capacity(degree + 1))
    }
}

/// Error that occurred during [`BSplineCurve`] construction because a wrong number of
/// control points was used.
#[derive(Debug, Clone)]
pub struct BSplineCurveError<N: Number, I: Infinitesimal<N>> {
    /// The B-spline that was used to construct the [`BSplineCurve`].
    pub bspline: BSpline<N>,
    /// The control points that was used to construct the [`BSplineCurve`].
    pub control_points: Vec<Jet<N, I>>,
}

/// A curve (more specific: a [spline]) that is constructed as a linear combination of a B-spline.
///
/// The underlying [`BSpline`] and control points can be either borrowed or owned, hence it has
/// either a static lifetime or the lifetime of its underlying data.
///
/// [spline]: https://en.wikipedia.org/wiki/Spline_(mathematics)
#[derive(Debug, Clone)]
pub struct BSplineCurve<N: Number, I: Infinitesimal<N>> {
    dim: Dim,
    bspline: BSpline<N>,
    control_points: Vec<Jet<N, I>>,
}

impl<N: Number, I: Infinitesimal<N>> BSplineCurve<N, I> {
    /// Tries to create a B-spline curve.
    ///
    /// The number of control points must match the number returned by
    /// [`BSpline::necessary_control_point_count`].
    pub fn new(
        dim: Dim,
        bspline: BSpline<N>,
        control_points: Vec<Jet<N, I>>,
    ) -> Result<BSplineCurve<N, I>, BSplineCurveError<N, I>> {
        if control_points.len() != bspline.necessary_control_point_count() {
            return Err(BSplineCurveError {
                bspline,
                control_points,
            });
        }

        Ok(BSplineCurve {
            bspline,
            dim,
            control_points,
        })
    }

    /// Returns the B-Spline of this curve.
    pub fn bspline(&self) -> &BSpline<N> {
        &self.bspline
    }

    /// Returns the control points of this curve.
    pub fn control_points(&self) -> &[Jet<N, I>] {
        &self.control_points
    }

    /// Evaluates the curve for the given x-value.
    pub fn value(&self, x: &N, buffer: &mut BSplineCurveBuffer<N, I>) -> Jet<N, I> {
        self.bspline.value(&self.control_points, buffer, x)
    }

    /// Evaluates the derivative of the curve for the given x-value.
    pub fn derivative(&self, x: &N, buffer: &mut BSplineCurveBuffer<N, I>) -> Jet<N, I> {
        self.bspline
            .derivative(self.dim, &self.control_points, buffer, x)
    }
}

fn min_knot_count(degree: usize) -> usize {
    degree + 1
}

fn necessary_control_point_count(degree: usize, knot_count: usize) -> usize {
    knot_count + degree.saturating_sub(1)
}

#[cfg(test)]
mod tests {
    use std::iter::repeat;

    use crate::bspline::{BSpline, BSplineCurve, BSplineCurveBuffer, BSplineErrorReason};
    use crate::{Dim, Jet, NoInfinitesimal};

    #[test]
    fn bspline_new_examples() {
        assert_eq!(
            BSpline::<f32>::new(0, vec![]).unwrap_err().reason,
            BSplineErrorReason::NotEnoughKnots,
        );

        assert_eq!(
            BSpline::<f32>::new(0, vec![1.0]).unwrap().padded_knots,
            [1.0]
        );

        assert_eq!(
            BSpline::<f32>::new(1, vec![1.0]).unwrap_err().reason,
            BSplineErrorReason::NotEnoughKnots,
        );

        assert_eq!(
            BSpline::<f32>::new(1, vec![1.0, 2.0]).unwrap().padded_knots,
            [1.0, 1.0, 2.0, 2.0]
        );

        assert_eq!(
            BSpline::<f32>::new(1, vec![1.0, 1.0]).unwrap_err().reason,
            BSplineErrorReason::NotStrictlyIncreasingKnots {
                wrong_knot_index: 1
            },
        );

        assert_eq!(
            BSpline::<f32>::new(1, vec![2.0, 1.0]).unwrap_err().reason,
            BSplineErrorReason::NotStrictlyIncreasingKnots {
                wrong_knot_index: 1
            },
        );

        assert_eq!(
            BSpline::<f32>::new(1, vec![f32::NAN, 1.0])
                .unwrap_err()
                .reason,
            BSplineErrorReason::NotStrictlyIncreasingKnots {
                wrong_knot_index: 1
            },
        );

        assert_eq!(
            BSpline::<f32>::new(1, vec![1.0, f32::NAN])
                .unwrap_err()
                .reason,
            BSplineErrorReason::NotStrictlyIncreasingKnots {
                wrong_knot_index: 1
            },
        );

        assert_eq!(
            BSpline::<f32>::new(1, vec![f32::NAN, f32::NAN])
                .unwrap_err()
                .reason,
            BSplineErrorReason::NotStrictlyIncreasingKnots {
                wrong_knot_index: 1
            },
        );

        assert_eq!(
            BSpline::<f32>::new(2, vec![1.0, 2.0, 3.0])
                .unwrap()
                .padded_knots,
            [1.0, 1.0, 1.0, 2.0, 3.0, 3.0, 3.0]
        );

        assert_eq!(
            BSpline::<f32>::new(3, vec![1.0, 2.0, 3.0, 4.0])
                .unwrap()
                .padded_knots,
            [1.0, 1.0, 1.0, 1.0, 2.0, 3.0, 4.0, 4.0, 4.0, 4.0]
        );
    }

    #[test]
    fn bspline_curve_new_example() {
        let dim = Dim(0);
        let control_points = |len: usize| {
            repeat(Jet::new(1.0 as f32, NoInfinitesimal))
                .take(len)
                .collect()
        };

        {
            let bspline = BSpline::<f32>::new(0, vec![1.0]).unwrap();

            assert!(BSplineCurve::new(dim, bspline.clone(), control_points(0)).is_err());
            assert!(BSplineCurve::new(dim, bspline.clone(), control_points(1)).is_ok());
            assert!(BSplineCurve::new(dim, bspline.clone(), control_points(2)).is_err());
        }

        {
            let bspline = BSpline::<f32>::new(0, vec![1.0, 2.0]).unwrap();

            assert!(BSplineCurve::new(dim, bspline.clone(), control_points(1)).is_err());
            assert!(BSplineCurve::new(dim, bspline.clone(), control_points(2)).is_ok());
            assert!(BSplineCurve::new(dim, bspline.clone(), control_points(3)).is_err());
        }

        {
            let bspline = BSpline::<f32>::new(1, vec![1.0, 2.0]).unwrap();

            assert!(BSplineCurve::new(dim, bspline.clone(), control_points(1)).is_err());
            assert!(BSplineCurve::new(dim, bspline.clone(), control_points(2)).is_ok());
            assert!(BSplineCurve::new(dim, bspline.clone(), control_points(3)).is_err());
        }

        {
            let bspline = BSpline::<f32>::new(1, vec![1.0, 2.0, 3.0]).unwrap();

            assert!(BSplineCurve::new(dim, bspline.clone(), control_points(2)).is_err());
            assert!(BSplineCurve::new(dim, bspline.clone(), control_points(3)).is_ok());
            assert!(BSplineCurve::new(dim, bspline.clone(), control_points(4)).is_err());
        }

        {
            let bspline = BSpline::<f32>::new(2, vec![1.0, 2.0, 3.0]).unwrap();

            assert!(BSplineCurve::new(dim, bspline.clone(), control_points(3)).is_err());
            assert!(BSplineCurve::new(dim, bspline.clone(), control_points(4)).is_ok());
            assert!(BSplineCurve::new(dim, bspline.clone(), control_points(5)).is_err());
        }

        {
            let bspline = BSpline::<f32>::new(3, vec![1.0, 2.0, 3.0, 4.0]).unwrap();

            assert!(BSplineCurve::new(dim, bspline.clone(), control_points(5)).is_err());
            assert!(BSplineCurve::new(dim, bspline.clone(), control_points(6)).is_ok());
            assert!(BSplineCurve::new(dim, bspline.clone(), control_points(7)).is_err());
        }
    }

    // Note: This B-spline curve is plotted in `bspline_examples.rs`
    #[test]
    fn bspline_example_0_a() {
        let degree = 0;
        let knots = vec![1.0];
        let dim = Dim(0);
        let control_points = [2.0].map(|p: f32| Jet::new(p, NoInfinitesimal)).to_vec();

        let bspline = BSpline::<f32>::new(degree, knots.clone()).unwrap();
        let curve = BSplineCurve::new(dim, bspline, control_points.clone()).unwrap();
        let mut buffer = BSplineCurveBuffer::new(degree);

        assert_eq!(curve.value(&0.5, &mut buffer).real, 2.0);
        assert_eq!(curve.value(&1.0, &mut buffer).real, 2.0);
        assert_eq!(curve.value(&1.5, &mut buffer).real, 2.0);

        assert_eq!(curve.derivative(&0.5, &mut buffer).real, 0.0);
        assert_eq!(curve.derivative(&1.0, &mut buffer).real, 0.0);
        assert_eq!(curve.derivative(&1.5, &mut buffer).real, 0.0);
    }

    // Note: This B-spline curve is plotted in `bspline_examples.rs`
    #[test]
    fn bspline_example_0_b() {
        let degree = 0;
        let knots = vec![1.0, 2.0, 3.0];
        let dim = Dim(0);
        let control_points = [2.0, 3.0, 1.0]
            .map(|p: f32| Jet::new(p, NoInfinitesimal))
            .to_vec();

        let bspline = BSpline::<f32>::new(degree, knots.clone()).unwrap();
        let curve = BSplineCurve::new(dim, bspline, control_points.clone()).unwrap();
        let mut buffer = BSplineCurveBuffer::new(degree);

        assert_eq!(curve.value(&0.5, &mut buffer).real, 2.0);
        assert_eq!(curve.value(&1.0, &mut buffer).real, 2.0);
        assert_eq!(curve.value(&1.5, &mut buffer).real, 2.0);
        assert_eq!(curve.value(&2.0, &mut buffer).real, 2.0);
        assert_eq!(curve.value(&2.5, &mut buffer).real, 3.0);
        assert_eq!(curve.value(&3.0, &mut buffer).real, 3.0);
        assert_eq!(curve.value(&3.5, &mut buffer).real, 1.0);

        assert_eq!(curve.derivative(&0.5, &mut buffer).real, 0.0);
        assert_eq!(curve.derivative(&1.0, &mut buffer).real, 0.0);
        assert_eq!(curve.derivative(&1.5, &mut buffer).real, 0.0);
        assert_eq!(curve.derivative(&2.0, &mut buffer).real, 0.0);
        assert_eq!(curve.derivative(&2.5, &mut buffer).real, 0.0);
        assert_eq!(curve.derivative(&3.0, &mut buffer).real, 0.0);
        assert_eq!(curve.derivative(&3.5, &mut buffer).real, 0.0);
    }

    // Note: This B-spline curve is plotted in `bspline_examples.rs`
    #[test]
    fn bspline_example_1_a() {
        let degree = 1;
        let knots = vec![1.0, 2.0];
        let dim = Dim(0);
        let control_points = [1.0, 3.0]
            .map(|p: f32| Jet::new(p, NoInfinitesimal))
            .to_vec();

        let bspline = BSpline::<f32>::new(degree, knots.clone()).unwrap();
        let curve = BSplineCurve::new(dim, bspline, control_points.clone()).unwrap();
        let mut buffer = BSplineCurveBuffer::new(degree);

        assert_eq!(curve.value(&0.5, &mut buffer).real, 0.0);
        assert_eq!(curve.value(&1.0, &mut buffer).real, 1.0);
        assert_eq!(curve.value(&1.5, &mut buffer).real, 2.0);
        assert_eq!(curve.value(&2.0, &mut buffer).real, 3.0);
        assert_eq!(curve.value(&2.5, &mut buffer).real, 4.0);

        assert_eq!(curve.derivative(&0.5, &mut buffer).real, 2.0);
        assert_eq!(curve.derivative(&1.0, &mut buffer).real, 2.0);
        assert_eq!(curve.derivative(&1.5, &mut buffer).real, 2.0);
        assert_eq!(curve.derivative(&2.0, &mut buffer).real, 2.0);
        assert_eq!(curve.derivative(&2.5, &mut buffer).real, 2.0);
    }

    // Note: This B-spline curve is plotted in `bspline_examples.rs`
    #[test]
    fn bspline_example_1_b() {
        let degree = 1;
        let knots = vec![1.0, 2.0, 3.0, 4.0];
        let dim = Dim(0);
        let control_points = [1.0, 3.0, 2.0, 5.0]
            .map(|p: f32| Jet::new(p, NoInfinitesimal))
            .to_vec();

        let bspline = BSpline::<f32>::new(degree, knots.clone()).unwrap();
        let curve = BSplineCurve::new(dim, bspline, control_points.clone()).unwrap();
        let mut buffer = BSplineCurveBuffer::new(degree);

        assert_eq!(curve.value(&0.5, &mut buffer).real, 0.0);
        assert_eq!(curve.value(&1.0, &mut buffer).real, 1.0);
        assert_eq!(curve.value(&1.5, &mut buffer).real, 2.0);
        assert_eq!(curve.value(&2.0, &mut buffer).real, 3.0);
        assert_eq!(curve.value(&2.5, &mut buffer).real, 2.5);
        assert_eq!(curve.value(&3.0, &mut buffer).real, 2.0);
        assert_eq!(curve.value(&3.5, &mut buffer).real, 3.5);
        assert_eq!(curve.value(&4.0, &mut buffer).real, 5.0);
        assert_eq!(curve.value(&4.5, &mut buffer).real, 6.5);

        assert_eq!(curve.derivative(&0.5, &mut buffer).real, 2.0);
        assert_eq!(curve.derivative(&1.0, &mut buffer).real, 2.0);
        assert_eq!(curve.derivative(&1.5, &mut buffer).real, 2.0);
        assert_eq!(curve.derivative(&2.0, &mut buffer).real, 2.0);
        assert_eq!(curve.derivative(&2.5, &mut buffer).real, -1.0);
        assert_eq!(curve.derivative(&3.0, &mut buffer).real, -1.0);
        assert_eq!(curve.derivative(&3.5, &mut buffer).real, 3.0);
        assert_eq!(curve.derivative(&4.0, &mut buffer).real, 3.0);
        assert_eq!(curve.derivative(&4.5, &mut buffer).real, 3.0);
    }

    // Note: This B-spline curve is plotted in `bspline_examples.rs`
    #[test]
    fn bspline_example_2_a() {
        let degree = 2;
        let knots = vec![1.0, 2.0, 3.0];
        let dim = Dim(0);
        let control_points = [1.0, 3.0, 2.0, 4.0]
            .map(|p: f32| Jet::new(p, NoInfinitesimal))
            .to_vec();

        let bspline = BSpline::<f32>::new(degree, knots.clone()).unwrap();
        let curve = BSplineCurve::new(dim, bspline, control_points.clone()).unwrap();
        let mut buffer = BSplineCurveBuffer::new(degree);

        assert_eq!(curve.value(&0.5, &mut buffer).real, -1.625);
        assert_eq!(curve.value(&1.0, &mut buffer).real, 1.0);
        assert_eq!(curve.value(&1.5, &mut buffer).real, 2.375);
        assert_eq!(curve.value(&2.0, &mut buffer).real, 2.5);
        assert_eq!(curve.value(&2.5, &mut buffer).real, 2.625);
        assert_eq!(curve.value(&3.0, &mut buffer).real, 4.0);
        assert_eq!(curve.value(&3.5, &mut buffer).real, 6.625);

        assert_eq!(curve.derivative(&0.5, &mut buffer).real, 6.5);
        assert_eq!(curve.derivative(&1.0, &mut buffer).real, 4.0);
        assert_eq!(curve.derivative(&1.5, &mut buffer).real, 1.5);
        assert_eq!(curve.derivative(&2.0, &mut buffer).real, -1.0);
        assert_eq!(curve.derivative(&2.5, &mut buffer).real, 1.5);
        assert_eq!(curve.derivative(&3.0, &mut buffer).real, 4.0);
        assert_eq!(curve.derivative(&3.5, &mut buffer).real, 6.5);
    }

    // Note: This B-spline curve is plotted in `bspline_examples.rs`
    #[test]
    fn bspline_example_2_b() {
        let degree = 2;
        let knots = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let dim = Dim(0);
        let control_points = [1.0, 3.0, 2.0, 4.0, 1.0, 5.0]
            .map(|p: f32| Jet::new(p, NoInfinitesimal))
            .to_vec();

        let bspline = BSpline::<f32>::new(degree, knots.clone()).unwrap();
        let curve = BSplineCurve::new(dim, bspline, control_points.clone()).unwrap();
        let mut buffer = BSplineCurveBuffer::new(degree);

        assert_eq!(curve.value(&0.5, &mut buffer).real, -1.625);
        assert_eq!(curve.value(&1.0, &mut buffer).real, 1.0);
        assert_eq!(curve.value(&1.5, &mut buffer).real, 2.375);
        assert_eq!(curve.value(&2.0, &mut buffer).real, 2.5);
        assert_eq!(curve.value(&2.5, &mut buffer).real, 2.375);
        assert_eq!(curve.value(&3.0, &mut buffer).real, 3.0);
        assert_eq!(curve.value(&3.5, &mut buffer).real, 3.375);
        assert_eq!(curve.value(&4.0, &mut buffer).real, 2.5);
        assert_eq!(curve.value(&4.5, &mut buffer).real, 2.375);
        assert_eq!(curve.value(&5.0, &mut buffer).real, 5.0);
        assert_eq!(curve.value(&5.5, &mut buffer).real, 10.375);

        assert_eq!(curve.derivative(&0.5, &mut buffer).real, 6.5);
        assert_eq!(curve.derivative(&1.0, &mut buffer).real, 4.0);
        assert_eq!(curve.derivative(&1.5, &mut buffer).real, 1.5);
        assert_eq!(curve.derivative(&2.0, &mut buffer).real, -1.0);
        assert_eq!(curve.derivative(&2.5, &mut buffer).real, 0.5);
        assert_eq!(curve.derivative(&3.0, &mut buffer).real, 2.0);
        assert_eq!(curve.derivative(&3.5, &mut buffer).real, -0.5);
        assert_eq!(curve.derivative(&4.0, &mut buffer).real, -3.0);
        assert_eq!(curve.derivative(&4.5, &mut buffer).real, 2.5);
        assert_eq!(curve.derivative(&5.0, &mut buffer).real, 8.0);
        assert_eq!(curve.derivative(&5.5, &mut buffer).real, 13.5);
    }

    // Note: This B-spline curve is plotted in `bspline_examples.rs`
    #[test]
    fn bspline_example_3_a() {
        let degree = 3;
        let knots = vec![1.0, 2.0, 3.0, 4.0];
        let dim = Dim(0);
        let control_points = [1.0, 3.0, 2.0, 4.0, 1.0, 5.0]
            .map(|p: f32| Jet::new(p, NoInfinitesimal))
            .to_vec();

        let bspline = BSpline::<f32>::new(degree, knots.clone()).unwrap();
        let curve = BSplineCurve::new(dim, bspline, control_points.clone()).unwrap();
        let mut buffer = BSplineCurveBuffer::new(degree);

        assert_eq!(curve.value(&0.5, &mut buffer).real, -4.260417);
        assert_eq!(curve.value(&1.0, &mut buffer).real, 1.0);
        assert_eq!(curve.value(&1.5, &mut buffer).real, 2.5104165);
        assert_eq!(curve.value(&2.0, &mut buffer).real, 2.5833333);
        assert_eq!(curve.value(&2.5, &mut buffer).real, 2.9375);
        assert_eq!(curve.value(&3.0, &mut buffer).real, 2.9166667);
        assert_eq!(curve.value(&3.5, &mut buffer).real, 2.3020833);
        assert_eq!(curve.value(&4.0, &mut buffer).real, 5.0);
        assert_eq!(curve.value(&4.5, &mut buffer).real, 15.947917);

        assert_eq!(curve.derivative(&0.5, &mut buffer).real, 15.8125);
        assert_eq!(curve.derivative(&1.0, &mut buffer).real, 6.0);
        assert_eq!(curve.derivative(&1.5, &mut buffer).real, 0.8125);
        assert_eq!(curve.derivative(&2.0, &mut buffer).real, 0.25);
        assert_eq!(curve.derivative(&2.5, &mut buffer).real, 0.75);
        assert_eq!(curve.derivative(&3.0, &mut buffer).real, -1.25);
        assert_eq!(curve.derivative(&3.5, &mut buffer).real, 0.4375);
        assert_eq!(curve.derivative(&4.0, &mut buffer).real, 12.0);
        assert_eq!(curve.derivative(&4.5, &mut buffer).real, 33.4375);
    }

    // Note: This B-spline curve is plotted in `bspline_examples.rs`
    #[test]
    fn bspline_example_3_b() {
        let degree = 3;
        let knots = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let dim = Dim(0);
        let control_points = [1.0, 3.0, 2.0, 4.0, 1.0, 5.0, 2.0, 3.0]
            .map(|p: f32| Jet::new(p, NoInfinitesimal))
            .to_vec();

        let bspline = BSpline::<f32>::new(degree, knots.clone()).unwrap();
        let curve = BSplineCurve::new(dim, bspline, control_points.clone()).unwrap();
        let mut buffer = BSplineCurveBuffer::new(degree);

        assert_eq!(curve.value(&0.5, &mut buffer).real, -4.260417);
        assert_eq!(curve.value(&1.0, &mut buffer).real, 1.0);
        assert_eq!(curve.value(&1.5, &mut buffer).real, 2.5104165);
        assert_eq!(curve.value(&2.0, &mut buffer).real, 2.5833333);
        assert_eq!(curve.value(&2.5, &mut buffer).real, 2.96875);
        assert_eq!(curve.value(&3.0, &mut buffer).real, 3.1666665);
        assert_eq!(curve.value(&3.5, &mut buffer).real, 2.5416665);
        assert_eq!(curve.value(&4.0, &mut buffer).real, 2.1666667);
        assert_eq!(curve.value(&4.5, &mut buffer).real, 2.96875);
        assert_eq!(curve.value(&5.0, &mut buffer).real, 3.5833335);
        assert_eq!(curve.value(&5.5, &mut buffer).real, 2.8854165);
        assert_eq!(curve.value(&6.0, &mut buffer).real, 3.0);
        assert_eq!(curve.value(&6.5, &mut buffer).real, 6.8645835);

        assert_eq!(curve.derivative(&0.5, &mut buffer).real, 15.8125);
        assert_eq!(curve.derivative(&1.0, &mut buffer).real, 6.0);
        assert_eq!(curve.derivative(&1.5, &mut buffer).real, 0.8125);
        assert_eq!(curve.derivative(&2.0, &mut buffer).real, 0.25);
        assert_eq!(curve.derivative(&2.5, &mut buffer).real, 0.9375);
        assert_eq!(curve.derivative(&3.0, &mut buffer).real, -0.5);
        assert_eq!(curve.derivative(&3.5, &mut buffer).real, -1.5);
        assert_eq!(curve.derivative(&4.0, &mut buffer).real, 0.5);
        assert_eq!(curve.derivative(&4.5, &mut buffer).real, 2.0625);
        assert_eq!(curve.derivative(&5.0, &mut buffer).real, -0.25);
        assert_eq!(curve.derivative(&5.5, &mut buffer).real, -1.5625);
        assert_eq!(curve.derivative(&6.0, &mut buffer).real, 3.0);
        assert_eq!(curve.derivative(&6.5, &mut buffer).real, 13.4375);
    }
}
