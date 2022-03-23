//! Provides an implementation for B-splines that uses jets as control points.
//!
//! Finding control points that fulfill the required properties is a non-trivial task.
//! One approach is to use a numeric solver to approximate the solution. With jets we
//! can calculate the derivatives with respect to the control points. This allows us to
//! use very efficient numeric solvers (e.g. Newton's method).

use std::{cmp::Ordering, fmt::Display, iter::repeat};

use crate::{Dim, Infinitesimal, Jet, NoInfinitesimal, Number};

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
        knot_count_without_padding(self.degree, &self.padded_knots)
    }

    /// The knots that were used to create the B-spline.
    pub fn knots(&self) -> &[N] {
        knots_without_padding(self.degree, &self.padded_knots)
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

impl<N: Number, I: Infinitesimal<N>> Default for BSplineCurveBuffer<N, I> {
    fn default() -> Self {
        Self::new(0)
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
    ) -> Result<Self, BSplineCurveError<N, I>> {
        if control_points.len() != bspline.necessary_control_point_count() {
            return Err(BSplineCurveError {
                bspline,
                control_points,
            });
        }

        Ok(Self {
            dim,
            bspline,
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
        calc_value(
            self.bspline.degree,
            &self.bspline.padded_knots,
            &self.control_points,
            x,
            buffer,
        )
    }

    /// Evaluates the curve's derivative of the given order and for the given x-value.
    pub fn derivative(
        &self,
        order: usize,
        x: &N,
        buffer: &mut BSplineCurveBuffer<N, I>,
    ) -> Jet<N, I> {
        calc_derivative(
            self.dim,
            self.bspline.degree,
            &self.bspline.padded_knots,
            &self.control_points,
            order,
            x,
            buffer,
        )
    }

    /// Returns a [`BSplineCurveRef`] that uses the underlying data of self.
    pub fn as_ref(&self) -> BSplineCurveRef<N, I> {
        BSplineCurveRef {
            dim: self.dim,
            bspline: &self.bspline,
            control_points: &self.control_points,
        }
    }
}

impl<N: Number> BSplineCurve<N, NoInfinitesimal> {
    /// Tries to create a B-spline curve using control points without an infinitesimal part.
    ///
    /// The number of control points must match the number returned by
    /// [`BSpline::necessary_control_point_count`].
    pub fn without_infinitesimal(
        bspline: BSpline<N>,
        real_control_points: Vec<N>,
    ) -> Result<Self, BSplineCurveError<N, NoInfinitesimal>> {
        let control_points = real_control_points
            .into_iter()
            .map(|p| Jet::new(p, NoInfinitesimal))
            .collect::<Vec<_>>();

        Self::new(Dim(0), bspline, control_points)
    }
}

/// Same as [`BSplineCurve`] but it borrows it underlying data.
#[derive(Debug, Clone)]
pub struct BSplineCurveRef<'a, N: Number, I: Infinitesimal<N>> {
    dim: Dim,
    bspline: &'a BSpline<N>,
    control_points: &'a [Jet<N, I>],
}

impl<'a, N: Number, I: Infinitesimal<N>> BSplineCurveRef<'a, N, I> {
    /// Tries to create a B-spline curve.
    ///
    /// The number of control points must match the number returned by
    /// [`BSpline::necessary_control_point_count`].
    pub fn new(dim: Dim, bspline: &'a BSpline<N>, control_points: &'a [Jet<N, I>]) -> Option<Self> {
        if control_points.len() != bspline.necessary_control_point_count() {
            return None;
        }

        Some(Self {
            bspline,
            dim,
            control_points,
        })
    }

    /// Returns the B-Spline of this curve.
    pub fn bspline(&self) -> &BSpline<N> {
        self.bspline
    }

    /// Returns the control points of this curve.
    pub fn control_points(&self) -> &[Jet<N, I>] {
        self.control_points
    }

    /// Evaluates the curve for the given x-value.
    pub fn value(&self, x: &N, buffer: &mut BSplineCurveBuffer<N, I>) -> Jet<N, I> {
        calc_value(
            self.bspline.degree,
            &self.bspline.padded_knots,
            self.control_points,
            x,
            buffer,
        )
    }

    /// Evaluates the curve's derivative of the given order and for the given x-value.
    pub fn derivative(
        &self,
        order: usize,
        x: &N,
        buffer: &mut BSplineCurveBuffer<N, I>,
    ) -> Jet<N, I> {
        calc_derivative(
            self.dim,
            self.bspline.degree,
            &self.bspline.padded_knots,
            self.control_points,
            order,
            x,
            buffer,
        )
    }
}

fn min_knot_count(degree: usize) -> usize {
    degree + 1
}

pub fn knot_count_without_padding<N: Number>(degree: usize, padded_knots: &[N]) -> usize {
    padded_knots.len() - 2 * degree
}

pub fn knots_without_padding<N: Number>(degree: usize, padded_knots: &[N]) -> &[N] {
    &padded_knots[degree..padded_knots.len() - degree]
}

fn necessary_control_point_count(degree: usize, knot_count: usize) -> usize {
    knot_count + degree.saturating_sub(1)
}

fn calc_value<N: Number, I: Infinitesimal<N>>(
    degree: usize,
    padded_knots: &[N],
    control_points: &[Jet<N, I>],
    x: &N,
    buffer: &mut BSplineCurveBuffer<N, I>,
) -> Jet<N, I> {
    // The algorithm is inspired by https://en.wikipedia.org/wiki/De_Boor%27s_algorithm

    // We assume that the number of control points correct because it was checked before
    debug_assert_eq!(
        control_points.len(),
        necessary_control_point_count(degree, knot_count_without_padding(degree, padded_knots)),
    );

    let interval_index = find_interval(degree, padded_knots, x);

    // Prepare the buffer with the relevant control points
    buffer.0.clear();
    buffer.0.extend(relevant_control_points(
        degree,
        interval_index,
        control_points,
    ));

    // Calculate the value
    de_boor_algorithm(degree, padded_knots, interval_index, x, buffer)
}

fn calc_derivative<N: Number, I: Infinitesimal<N>>(
    dim: Dim,
    degree: usize,
    padded_knots: &[N],
    control_points: &[Jet<N, I>],
    order: usize,
    x: &N,
    buffer: &mut BSplineCurveBuffer<N, I>,
) -> Jet<N, I> {
    if order == 0 {
        // We interpret the derivative of order 0 as the original curve
        return calc_value(degree, padded_knots, control_points, x, buffer);
    }

    // We assume that the number of control points correct because it was checked before
    debug_assert_eq!(
        control_points.len(),
        necessary_control_point_count(degree, knot_count_without_padding(degree, padded_knots)),
    );

    if order > degree {
        // If the order is greater than the degree, we knot that derivative is always zero
        return Jet::new(N::zero(), I::zeros(dim));
    }

    let derivative_degree = degree - order;
    let interval_index = find_interval(degree, padded_knots, x);

    // Prepare the buffer with the relevant control points
    buffer.0.clear();
    buffer.0.extend(relevant_control_points(
        degree,
        interval_index,
        control_points,
    ));

    // Calculate the control points for the derivative
    calc_derivative_control_points(
        degree,
        derivative_degree,
        padded_knots,
        interval_index,
        buffer,
    );

    // Calculate the derivative
    de_boor_algorithm(derivative_degree, padded_knots, interval_index, x, buffer)
}

fn find_interval<N: Number>(degree: usize, padded_knots: &[N], x: &N) -> usize {
    let knots = knots_without_padding(degree, padded_knots);

    // Calculate the knots between the intervals
    let interval_knots = if degree == 0 {
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

    interval_knot_index + degree
}

fn relevant_control_points<N: Number, I: Infinitesimal<N>>(
    degree: usize,
    interval_index: usize,
    control_points: &[Jet<N, I>],
) -> impl Iterator<Item = Jet<N, I>> + '_ {
    control_points[interval_index - degree..=interval_index]
        .iter()
        .cloned()
}

// Implements the recursive part of the De Boor's algorithm.
//
// Expects that the buffer contains `degree + 1` control points.
//
// Inspired by https://en.wikipedia.org/wiki/De_Boor%27s_algorithm
fn de_boor_algorithm<N: Number, I: Infinitesimal<N>>(
    degree: usize,
    padded_knots: &[N],
    interval_index: usize,
    x: &N,
    buffer: &mut BSplineCurveBuffer<N, I>,
) -> Jet<N, I> {
    for step in 1..=degree {
        for i in (step..=degree).rev() {
            let left_knot = &padded_knots[interval_index + i - degree];
            let right_knot = &padded_knots[interval_index + 1 + i - step];
            let knot_diff = right_knot.clone() - left_knot.clone();

            // `knot_diff` is non-zero because we assume that the knots are strictly increasing
            // and at most one padding knot is used for the difference
            debug_assert!(!knot_diff.is_zero());

            let alpha = (x.clone() - left_knot.clone()) / knot_diff;
            let beta = N::one() - alpha.clone();

            let left_temp = buffer.0[i - 1].clone();
            let right_temp = &mut buffer.0[i];
            *right_temp = left_temp * beta + right_temp.clone() * alpha;
        }
    }

    buffer.0.pop().unwrap()
}

// Calculates the control points for the curve's derivative with the given degree.
//
// Expects that the buffer contains `degree + 1` control points of the B-spline curve.
// After this functions returns the buffer contains `derivative_degree + 1` control points of
// the curve's derivative.
//
// Inspired by https://stackoverflow.com/questions/57507696/b-spline-derivative-using-de-boors-algorithm
fn calc_derivative_control_points<N: Number, I: Infinitesimal<N>>(
    degree: usize,
    derivative_degree: usize,
    padded_knots: &[N],
    interval_index: usize,
    buffer: &mut BSplineCurveBuffer<N, I>,
) {
    for step in (derivative_degree..degree).rev() {
        for i in 0..=step {
            let left_knot = &padded_knots[i + interval_index - step];
            let right_knot = &padded_knots[i + interval_index + 1];
            let knot_diff = right_knot.clone() - left_knot.clone();

            // `knot_diff` is non-zero because we assume that the knots are strictly increasing
            // and at most one padding knot is used for the difference
            debug_assert!(!knot_diff.is_zero());

            let alpha = N::from_integer((step + 1) as i32) / knot_diff;

            let right_temp = buffer.0[i + 1].clone();
            let left_temp = &mut buffer.0[i];
            *left_temp = (right_temp - left_temp.clone()) * alpha;
        }
    }

    // Remove unnecessary entries from buffer
    buffer.0.drain(derivative_degree + 1..);
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

    fn create_bspline_example(
        degree: usize,
        knots: Vec<f32>,
        control_points: Vec<f32>,
    ) -> BSplineCurve<f32, NoInfinitesimal> {
        let dim = Dim(0);
        let control_points = control_points
            .into_iter()
            .map(|p| Jet::new(p, NoInfinitesimal))
            .collect();
        let bspline = BSpline::<f32>::new(degree, knots).unwrap();
        BSplineCurve::new(dim, bspline, control_points).unwrap()
    }

    fn check_bspline_example(
        curve: BSplineCurve<f32, NoInfinitesimal>,
        expected_results: Vec<(Option<usize>, f32, f32)>,
    ) {
        let mut violations = Vec::new();
        let mut buffer = BSplineCurveBuffer::new(curve.bspline().degree);

        for (order, x, expected_y) in expected_results {
            let actual_y = match order {
                None => curve.value(&x, &mut buffer).real,
                Some(order) => curve.derivative(order, &x, &mut buffer).real,
            };

            if expected_y != actual_y {
                violations.push((order, x, expected_y, actual_y))
            }
        }

        let no_violation = violations.is_empty();

        let violation_message = move || {
            let mut message = String::new();

            message += "B-spline produces unexpected results:\n";

            for (order, x, expected_y, actual_y) in violations {
                let prefix = match order {
                    None => String::from("f"),
                    Some(order) => format!("f^{order}"),
                };
                message +=
                    &format!("\t{prefix}({x}) is {actual_y}, but {expected_y} was expected\n");
            }

            message
        };

        assert!(no_violation, "{}", violation_message());
    }

    // Note: This B-spline curve is plotted in `bspline_examples.rs`
    #[test]
    fn bspline_example_0_a() {
        let curve = create_bspline_example(0, vec![1.0], vec![2.0]);

        check_bspline_example(
            curve,
            vec![
                (None, 0.5, 2.0),
                (None, 1.0, 2.0),
                (None, 1.5, 2.0),
                (Some(0), 0.5, 2.0),
                (Some(0), 1.0, 2.0),
                (Some(0), 1.5, 2.0),
                (Some(1), 0.5, 0.0),
                (Some(1), 1.0, 0.0),
                (Some(1), 1.5, 0.0),
                (Some(2), 0.5, 0.0),
                (Some(2), 1.0, 0.0),
                (Some(2), 1.5, 0.0),
            ],
        );
    }

    // Note: This B-spline curve is plotted in `bspline_examples.rs`
    #[test]
    fn bspline_example_0_b() {
        let curve = create_bspline_example(0, vec![1.0, 2.0, 3.0], vec![2.0, 3.0, 1.0]);

        check_bspline_example(
            curve,
            vec![
                (None, 0.5, 2.0),
                (None, 1.0, 2.0),
                (None, 1.5, 2.0),
                (None, 2.0, 2.0),
                (None, 2.5, 3.0),
                (None, 3.0, 3.0),
                (None, 3.5, 1.0),
                (Some(0), 0.5, 2.0),
                (Some(0), 1.0, 2.0),
                (Some(0), 1.5, 2.0),
                (Some(0), 2.0, 2.0),
                (Some(0), 2.5, 3.0),
                (Some(0), 3.0, 3.0),
                (Some(0), 3.5, 1.0),
                (Some(1), 0.5, 0.0),
                (Some(1), 1.0, 0.0),
                (Some(1), 1.5, 0.0),
                (Some(1), 2.0, 0.0),
                (Some(1), 2.5, 0.0),
                (Some(1), 3.0, 0.0),
                (Some(1), 3.5, 0.0),
                (Some(2), 0.5, 0.0),
                (Some(2), 1.0, 0.0),
                (Some(2), 1.5, 0.0),
                (Some(2), 2.0, 0.0),
                (Some(2), 2.5, 0.0),
                (Some(2), 3.0, 0.0),
                (Some(2), 3.5, 0.0),
            ],
        );
    }

    // Note: This B-spline curve is plotted in `bspline_examples.rs`
    #[test]
    fn bspline_example_1_a() {
        let curve = create_bspline_example(1, vec![1.0, 2.0], vec![1.0, 3.0]);

        check_bspline_example(
            curve,
            vec![
                (None, 0.5, 0.0),
                (None, 1.0, 1.0),
                (None, 1.5, 2.0),
                (None, 2.0, 3.0),
                (None, 2.5, 4.0),
                (Some(0), 0.5, 0.0),
                (Some(0), 1.0, 1.0),
                (Some(0), 1.5, 2.0),
                (Some(0), 2.0, 3.0),
                (Some(0), 2.5, 4.0),
                (Some(1), 0.5, 2.0),
                (Some(1), 1.0, 2.0),
                (Some(1), 1.5, 2.0),
                (Some(1), 2.0, 2.0),
                (Some(1), 2.5, 2.0),
                (Some(2), 0.5, 0.0),
                (Some(2), 1.0, 0.0),
                (Some(2), 1.5, 0.0),
                (Some(2), 2.0, 0.0),
                (Some(2), 2.5, 0.0),
                (Some(3), 0.5, 0.0),
                (Some(3), 1.0, 0.0),
                (Some(3), 1.5, 0.0),
                (Some(3), 2.0, 0.0),
                (Some(3), 2.5, 0.0),
            ],
        );
    }

    // Note: This B-spline curve is plotted in `bspline_examples.rs`
    #[test]
    fn bspline_example_1_b() {
        let curve = create_bspline_example(1, vec![1.0, 2.0, 3.0, 4.0], vec![1.0, 3.0, 2.0, 5.0]);

        check_bspline_example(
            curve,
            vec![
                (None, 0.5, 0.0),
                (None, 1.0, 1.0),
                (None, 1.5, 2.0),
                (None, 2.0, 3.0),
                (None, 2.5, 2.5),
                (None, 3.0, 2.0),
                (None, 3.5, 3.5),
                (None, 4.0, 5.0),
                (None, 4.5, 6.5),
                (Some(0), 0.5, 0.0),
                (Some(0), 1.0, 1.0),
                (Some(0), 1.5, 2.0),
                (Some(0), 2.0, 3.0),
                (Some(0), 2.5, 2.5),
                (Some(0), 3.0, 2.0),
                (Some(0), 3.5, 3.5),
                (Some(0), 4.0, 5.0),
                (Some(0), 4.5, 6.5),
                (Some(1), 0.5, 2.0),
                (Some(1), 1.0, 2.0),
                (Some(1), 1.5, 2.0),
                (Some(1), 2.0, 2.0),
                (Some(1), 2.5, -1.0),
                (Some(1), 3.0, -1.0),
                (Some(1), 3.5, 3.0),
                (Some(1), 4.0, 3.0),
                (Some(1), 4.5, 3.0),
                (Some(2), 0.5, 0.0),
                (Some(2), 1.0, 0.0),
                (Some(2), 1.5, 0.0),
                (Some(2), 2.0, 0.0),
                (Some(2), 2.5, 0.0),
                (Some(2), 3.0, 0.0),
                (Some(2), 3.5, 0.0),
                (Some(2), 4.0, 0.0),
                (Some(2), 4.5, 0.0),
                (Some(3), 0.5, 0.0),
                (Some(3), 1.0, 0.0),
                (Some(3), 1.5, 0.0),
                (Some(3), 2.0, 0.0),
                (Some(3), 2.5, 0.0),
                (Some(3), 3.0, 0.0),
                (Some(3), 3.5, 0.0),
                (Some(3), 4.0, 0.0),
                (Some(3), 4.5, 0.0),
            ],
        );
    }

    // Note: This B-spline curve is plotted in `bspline_examples.rs`
    #[test]
    fn bspline_example_2_a() {
        let curve = create_bspline_example(2, vec![1.0, 2.0, 3.0], vec![1.0, 3.0, 2.0, 4.0]);

        check_bspline_example(
            curve,
            vec![
                (None, 0.5, -1.625),
                (None, 1.0, 1.0),
                (None, 1.5, 2.375),
                (None, 2.0, 2.5),
                (None, 2.5, 2.625),
                (None, 3.0, 4.0),
                (None, 3.5, 6.625),
                (Some(0), 0.5, -1.625),
                (Some(0), 1.0, 1.0),
                (Some(0), 1.5, 2.375),
                (Some(0), 2.0, 2.5),
                (Some(0), 2.5, 2.625),
                (Some(0), 3.0, 4.0),
                (Some(0), 3.5, 6.625),
                (Some(1), 0.5, 6.5),
                (Some(1), 1.0, 4.0),
                (Some(1), 1.5, 1.5),
                (Some(1), 2.0, -1.0),
                (Some(1), 2.5, 1.5),
                (Some(1), 3.0, 4.0),
                (Some(1), 3.5, 6.5),
                (Some(2), 0.5, -5.0),
                (Some(2), 1.0, -5.0),
                (Some(2), 1.5, -5.0),
                (Some(2), 2.0, -5.0),
                (Some(2), 2.5, 5.0),
                (Some(2), 3.0, 5.0),
                (Some(2), 3.5, 5.0),
                (Some(3), 0.5, 0.0),
                (Some(3), 1.0, 0.0),
                (Some(3), 1.5, 0.0),
                (Some(3), 2.0, 0.0),
                (Some(3), 2.5, 0.0),
                (Some(3), 3.0, 0.0),
                (Some(3), 3.5, 0.0),
                (Some(4), 0.5, 0.0),
                (Some(4), 1.0, 0.0),
                (Some(4), 1.5, 0.0),
                (Some(4), 2.0, 0.0),
                (Some(4), 2.5, 0.0),
                (Some(4), 3.0, 0.0),
                (Some(4), 3.5, 0.0),
            ],
        );
    }

    // Note: This B-spline curve is plotted in `bspline_examples.rs`
    #[test]
    fn bspline_example_2_b() {
        let curve = create_bspline_example(
            2,
            vec![1.0, 2.0, 3.0, 4.0, 5.0],
            vec![1.0, 3.0, 2.0, 4.0, 1.0, 5.0],
        );

        check_bspline_example(
            curve,
            vec![
                (None, 0.5, -1.625),
                (None, 1.0, 1.0),
                (None, 1.5, 2.375),
                (None, 2.0, 2.5),
                (None, 2.5, 2.375),
                (None, 3.0, 3.0),
                (None, 3.5, 3.375),
                (None, 4.0, 2.5),
                (None, 4.5, 2.375),
                (None, 5.0, 5.0),
                (None, 5.5, 10.375),
                (Some(0), 0.5, -1.625),
                (Some(0), 1.0, 1.0),
                (Some(0), 1.5, 2.375),
                (Some(0), 2.0, 2.5),
                (Some(0), 2.5, 2.375),
                (Some(0), 3.0, 3.0),
                (Some(0), 3.5, 3.375),
                (Some(0), 4.0, 2.5),
                (Some(0), 4.5, 2.375),
                (Some(0), 5.0, 5.0),
                (Some(0), 5.5, 10.375),
                (Some(1), 0.5, 6.5),
                (Some(1), 1.0, 4.0),
                (Some(1), 1.5, 1.5),
                (Some(1), 2.0, -1.0),
                (Some(1), 2.5, 0.5),
                (Some(1), 3.0, 2.0),
                (Some(1), 3.5, -0.5),
                (Some(1), 4.0, -3.0),
                (Some(1), 4.5, 2.5),
                (Some(1), 5.0, 8.0),
                (Some(1), 5.5, 13.5),
                (Some(2), 0.5, -5.0),
                (Some(2), 1.0, -5.0),
                (Some(2), 1.5, -5.0),
                (Some(2), 2.0, -5.0),
                (Some(2), 2.5, 3.0),
                (Some(2), 3.0, 3.0),
                (Some(2), 3.5, -5.0),
                (Some(2), 4.0, -5.0),
                (Some(2), 4.5, 11.0),
                (Some(2), 5.0, 11.0),
                (Some(2), 5.5, 11.0),
                (Some(3), 0.5, 0.0),
                (Some(3), 1.0, 0.0),
                (Some(3), 1.5, 0.0),
                (Some(3), 2.0, 0.0),
                (Some(3), 2.5, 0.0),
                (Some(3), 3.0, 0.0),
                (Some(3), 3.5, 0.0),
                (Some(3), 4.0, 0.0),
                (Some(3), 4.5, 0.0),
                (Some(3), 5.0, 0.0),
                (Some(3), 5.5, 0.0),
                (Some(4), 0.5, 0.0),
                (Some(4), 1.0, 0.0),
                (Some(4), 1.5, 0.0),
                (Some(4), 2.0, 0.0),
                (Some(4), 2.5, 0.0),
                (Some(4), 3.0, 0.0),
                (Some(4), 3.5, 0.0),
                (Some(4), 4.0, 0.0),
                (Some(4), 4.5, 0.0),
                (Some(4), 5.0, 0.0),
                (Some(4), 5.5, 0.0),
            ],
        );
    }

    // Note: This B-spline curve is plotted in `bspline_examples.rs`
    #[test]
    fn bspline_example_3_a() {
        let curve = create_bspline_example(
            3,
            vec![1.0, 2.0, 3.0, 4.0],
            vec![1.0, 3.0, 2.0, 4.0, 1.0, 5.0],
        );

        check_bspline_example(
            curve,
            vec![
                (None, 0.5, -4.260417),
                (None, 1.0, 1.0),
                (None, 1.5, 2.5104165),
                (None, 2.0, 2.5833333),
                (None, 2.5, 2.9375),
                (None, 3.0, 2.9166667),
                (None, 3.5, 2.3020833),
                (None, 4.0, 5.0),
                (None, 4.5, 15.947917),
                (Some(0), 0.5, -4.260417),
                (Some(0), 1.0, 1.0),
                (Some(0), 1.5, 2.5104165),
                (Some(0), 2.0, 2.5833333),
                (Some(0), 2.5, 2.9375),
                (Some(0), 3.0, 2.9166667),
                (Some(0), 3.5, 2.3020833),
                (Some(0), 4.0, 5.0),
                (Some(0), 4.5, 15.947917),
                (Some(1), 0.5, 15.8125),
                (Some(1), 1.0, 6.0),
                (Some(1), 1.5, 0.8125),
                (Some(1), 2.0, 0.25),
                (Some(1), 2.5, 0.75),
                (Some(1), 3.0, -1.25),
                (Some(1), 3.5, 0.4375),
                (Some(1), 4.0, 12.0),
                (Some(1), 4.5, 33.4375),
                (Some(2), 0.5, -24.25),
                (Some(2), 1.0, -15.0),
                (Some(2), 1.5, -5.75),
                (Some(2), 2.0, 3.5),
                (Some(2), 2.5, -1.5),
                (Some(2), 3.0, -6.5),
                (Some(2), 3.5, 13.25),
                (Some(2), 4.0, 33.0),
                (Some(2), 4.5, 52.75),
                (Some(3), 0.5, 18.5),
                (Some(3), 1.0, 18.50),
                (Some(3), 1.5, 18.5),
                (Some(3), 2.0, 18.5),
                (Some(3), 2.5, -10.0),
                (Some(3), 3.0, -10.0),
                (Some(3), 3.5, 39.5),
                (Some(3), 4.0, 39.5),
                (Some(3), 4.5, 39.5),
                (Some(4), 0.5, 0.0),
                (Some(4), 1.0, 0.0),
                (Some(4), 1.5, 0.0),
                (Some(4), 2.0, 0.0),
                (Some(4), 2.5, 0.0),
                (Some(4), 3.0, 0.0),
                (Some(4), 3.5, 0.0),
                (Some(4), 4.0, 0.0),
                (Some(4), 4.5, 0.0),
                (Some(5), 0.5, 0.0),
                (Some(5), 1.0, 0.0),
                (Some(5), 1.5, 0.0),
                (Some(5), 2.0, 0.0),
                (Some(5), 2.5, 0.0),
                (Some(5), 3.0, 0.0),
                (Some(5), 3.5, 0.0),
                (Some(5), 4.0, 0.0),
                (Some(5), 4.5, 0.0),
            ],
        );
    }

    // Note: This B-spline curve is plotted in `bspline_examples.rs`
    #[test]
    fn bspline_example_3_b() {
        let curve = create_bspline_example(
            3,
            vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            vec![1.0, 3.0, 2.0, 4.0, 1.0, 5.0, 2.0, 3.0],
        );

        check_bspline_example(
            curve,
            vec![
                (None, 0.5, -4.260417),
                (None, 1.0, 1.0),
                (None, 1.5, 2.5104165),
                (None, 2.0, 2.5833333),
                (None, 2.5, 2.96875),
                (None, 3.0, 3.1666665),
                (None, 3.5, 2.5416665),
                (None, 4.0, 2.1666667),
                (None, 4.5, 2.96875),
                (None, 5.0, 3.5833335),
                (None, 5.5, 2.8854165),
                (None, 6.0, 3.0),
                (None, 6.5, 6.8645835),
                (Some(0), 0.5, -4.260417),
                (Some(0), 1.0, 1.0),
                (Some(0), 1.5, 2.5104165),
                (Some(0), 2.0, 2.5833333),
                (Some(0), 2.5, 2.96875),
                (Some(0), 3.0, 3.1666665),
                (Some(0), 3.5, 2.5416665),
                (Some(0), 4.0, 2.1666667),
                (Some(0), 4.5, 2.96875),
                (Some(0), 5.0, 3.5833335),
                (Some(0), 5.5, 2.8854165),
                (Some(0), 6.0, 3.0),
                (Some(0), 6.5, 6.8645835),
                (Some(1), 0.5, 15.8125),
                (Some(1), 1.0, 6.0),
                (Some(1), 1.5, 0.8125),
                (Some(1), 2.0, 0.25),
                (Some(1), 2.5, 0.9375),
                (Some(1), 3.0, -0.5),
                (Some(1), 3.5, -1.5),
                (Some(1), 4.0, 0.5),
                (Some(1), 4.5, 2.0625),
                (Some(1), 5.0, -0.25),
                (Some(1), 5.5, -1.5625),
                (Some(1), 6.0, 3.0),
                (Some(1), 6.5, 13.4375),
                (Some(2), 0.5, -24.25),
                (Some(2), 1.0, -15.0),
                (Some(2), 1.5, -5.75),
                (Some(2), 2.0, 3.5),
                (Some(2), 2.5, -0.75),
                (Some(2), 3.0, -5.0),
                (Some(2), 3.5, 1.0),
                (Some(2), 4.0, 7.0),
                (Some(2), 4.5, -0.75),
                (Some(2), 5.0, -8.5),
                (Some(2), 5.5, 3.25),
                (Some(2), 6.0, 15.0),
                (Some(2), 6.5, 26.75),
                (Some(3), 0.5, 18.5),
                (Some(3), 1.0, 18.5),
                (Some(3), 1.5, 18.5),
                (Some(3), 2.0, 18.5),
                (Some(3), 2.5, -8.5),
                (Some(3), 3.0, -8.5),
                (Some(3), 3.5, 12.0),
                (Some(3), 4.0, 12.0),
                (Some(3), 4.5, -15.5),
                (Some(3), 5.0, -15.5),
                (Some(3), 5.5, 23.5),
                (Some(3), 6.0, 23.5),
                (Some(3), 6.5, 23.5),
                (Some(4), 0.5, 0.0),
                (Some(4), 1.0, 0.0),
                (Some(4), 1.5, 0.0),
                (Some(4), 2.0, 0.0),
                (Some(4), 2.5, 0.0),
                (Some(4), 3.0, 0.0),
                (Some(4), 3.5, 0.0),
                (Some(4), 4.0, 0.0),
                (Some(4), 4.5, 0.0),
                (Some(4), 5.0, 0.0),
                (Some(4), 5.5, 0.0),
                (Some(4), 6.0, 0.0),
                (Some(4), 6.5, 0.0),
                (Some(5), 0.5, 0.0),
                (Some(5), 1.0, 0.0),
                (Some(5), 1.5, 0.0),
                (Some(5), 2.0, 0.0),
                (Some(5), 2.5, 0.0),
                (Some(5), 3.0, 0.0),
                (Some(5), 3.5, 0.0),
                (Some(5), 4.0, 0.0),
                (Some(5), 4.5, 0.0),
                (Some(5), 5.0, 0.0),
                (Some(5), 5.5, 0.0),
                (Some(5), 6.0, 0.0),
                (Some(5), 6.5, 0.0),
            ],
        );
    }
}
