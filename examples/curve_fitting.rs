// This example demonstrates how jets can be used for curve fitting.

use std::iter::repeat;
use std::path::Path;
use std::{fs::File, io::Write};

use poloto::prelude::*;
use wtfisjet::bspline::{BSpline, BSplineCurve, BSplineCurveBuffer, BSplineCurveRef};
use wtfisjet::newton::NewtonsMethod;
use wtfisjet::{DenseInfinitesimal, NoInfinitesimal};

fn main() {
    // The degree of the B-Spline
    let degree = 3;
    // The knots of the B-Spline
    let knots = vec![
        1.0, 2.0, 4.0, 6.0, 9.0, 10.0, 13.0, 15.0, 17.0, 19.0, 25.0, 26.0, 28.0, 29.0,
    ];
    // The points that the B-Spline curve must match
    let curve_constraints = vec![
        (1.0, 4.0),
        (2.0, -4.0),
        (3.0, 10.0),
        (7.0, 15.0),
        (9.0, -6.0),
        (11.0, 5.0),
        (14.0, -10.0),
        (16.0, 5.0),
        (17.0, 3.0),
        (22.0, 10.0),
        (23.0, 0.0),
        (25.0, -12.0),
        (28.0, 19.0),
        (29.0, 9.0),
    ];
    // The points that the derivative of the B-Spline curve must match
    let derivative_constraints = vec![(1.0, -5.0), (29.0, -10.0)];

    // Run the curve fitting
    let curve = fit_curve(
        degree,
        knots.clone(),
        curve_constraints.clone(),
        derivative_constraints.clone(),
    );

    // Draw the result
    plot_curve(
        "curve fitting",
        "curve_fitting.svg",
        knots,
        curve_constraints,
        derivative_constraints,
        curve,
    );
}

fn fit_curve(
    degree: usize,
    knots: Vec<f64>,
    curve_constraints: Vec<(f64, f64)>,
    derivative_constraints: Vec<(f64, f64)>,
) -> BSplineCurve<f64, NoInfinitesimal> {
    let bspline = BSpline::new(degree, knots).expect("Failed to create the B-spline");
    let mut buffer = BSplineCurveBuffer::new(degree);

    let newtons_method = NewtonsMethod::<f64> {
        max_steps: 10,
        threshold: 1e-10,
    };

    // The control points of the B-spline are the variables of the root-finder.
    // The error (difference between expected and actual points) is the root that we are searching.
    let result = newtons_method.find_root::<DenseInfinitesimal<f64>>(
        repeat(1.0)
            .take(bspline.necessary_control_point_count())
            .collect(),
        &mut |dim, control_points, errors| {
            let curve = BSplineCurveRef::new(dim, &bspline, control_points)
                .expect("Failed to create the temporary B-spline curve");

            for (x, y) in curve_constraints.iter() {
                let current_y = curve.value(x, &mut buffer);
                let pillar_error = current_y - *y;
                errors.push(pillar_error);
            }

            for (x, y) in derivative_constraints.iter() {
                let current_y = curve.derivative(1, x, &mut buffer);
                let pillar_error = current_y - *y;
                errors.push(pillar_error);
            }
        },
    );

    println!("result = {result:?}");

    // Construct a curve using the latest guess of the root-finder
    BSplineCurve::without_infinitesimal(bspline, result.latest_guess)
        .expect("Failed to create the final B-spline curve")
}

fn plot_curve(
    title_name: &str,
    file_name: &str,
    knots: Vec<f64>,
    curve_constraints: Vec<(f64, f64)>,
    derivative_constraints: Vec<(f64, f64)>,
    curve: BSplineCurve<f64, NoInfinitesimal>,
) {
    // Create the plot
    let knots_points = knots
        .into_iter()
        .map(|k| [k, 0.0 as f64])
        .collect::<Vec<_>>();
    let knots_plot = poloto::build::scatter("knot", knots_points);

    let curve_constraints_plot = poloto::build::scatter("f constraint", curve_constraints);

    let derivative_constraints_plot =
        poloto::build::scatter("f' constraint", derivative_constraints);

    let mut buffer = BSplineCurveBuffer::new(curve.bspline().degree());
    let range = {
        let knots = curve.bspline().knots();
        let min = knots.first().unwrap() - 1.0;
        let max = knots.last().unwrap() + 1.0;
        let delta = max - min;
        let points = 1000;
        (0..points).map(move |x| min + (x as f64 * delta) / points as f64)
    };

    let curve_points = range
        .clone()
        .map(|x| [x as f64, curve.value(&x, &mut buffer).real as f64])
        .collect::<Vec<_>>();
    let curve_plot = poloto::build::line("f", curve_points);

    let derivative_points = range
        .clone()
        .map(|x| [x as f64, curve.derivative(1, &x, &mut buffer).real as f64])
        .collect::<Vec<_>>();
    let derivative_plot = poloto::build::line("f'", derivative_points);

    let plots = plots!(
        knots_plot,
        derivative_constraints_plot,
        curve_constraints_plot,
        derivative_plot,
        curve_plot
    )
    .markers([], [0.0]);

    let plotter = simple_fmt!(plots, title_name, "x", "y");

    // Write the plot to a file
    let folder = Path::new("plots");
    let file = folder.join(file_name);

    std::fs::create_dir_all(folder).expect("Failed to create the folder for the plots");

    let mut file = File::create(file).expect("Failed to create the file for plot");

    write!(file, "{}", poloto::disp(|f| plotter.simple_theme(f)))
        .expect("Failed to create write the plot");
}
