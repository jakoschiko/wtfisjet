// Plots the example B-splines from the unit tests.

use std::path::Path;
use std::{fs::File, io::Write};

use poloto::prelude::*;
use wtfisjet::bspline::{BSpline, BSplineCurve, BSplineCurveBuffer};
use wtfisjet::{Dim, Jet, NoInfinitesimal};

fn main() {
    // Create the B-spline curves
    let curve_0_a = bspline_example_0_a();
    let curve_0_b = bspline_example_0_b();
    let curve_1_a = bspline_example_1_a();
    let curve_1_b = bspline_example_1_b();
    let curve_2_a = bspline_example_2_a();
    let curve_2_b = bspline_example_2_b();
    let curve_3_a = bspline_example_3_a();
    let curve_3_b = bspline_example_3_b();

    // Plot the curves
    plot_curve("bspline example 0 a", "bspline_example_0_a.svg", curve_0_a);
    plot_curve("bspline example 0 b", "bspline_example_0_b.svg", curve_0_b);
    plot_curve("bspline example 1 a", "bspline_example_1_a.svg", curve_1_a);
    plot_curve("bspline example 1 b", "bspline_example_1_b.svg", curve_1_b);
    plot_curve("bspline example 2 a", "bspline_example_2_a.svg", curve_2_a);
    plot_curve("bspline example 2 b", "bspline_example_2_b.svg", curve_2_b);
    plot_curve("bspline example 3 a", "bspline_example_3_a.svg", curve_3_a);
    plot_curve("bspline example 3 b", "bspline_example_3_b.svg", curve_3_b);
}

fn plot_curve(title_name: &str, file_name: &str, curve: BSplineCurve<f32, NoInfinitesimal>) {
    // Create the plot
    let knots_points = curve
        .bspline()
        .knots()
        .iter()
        .map(|&k| [k as f64, 0.0 as f64])
        .collect::<Vec<_>>();
    let knots_plot = poloto::build::scatter("knot", knots_points);

    let mut buffer = BSplineCurveBuffer::new(curve.bspline().degree());
    let range = {
        let knots = curve.bspline().knots();
        let min = knots.first().unwrap() - 1.0;
        let max = knots.last().unwrap() + 1.0;
        let delta = max - min;
        let points = 1000;
        (0..points).map(move |x| min + (x as f32 * delta) / points as f32)
    };

    let curve_points = range
        .clone()
        .map(|x| [x as f64, curve.value(&x, &mut buffer).real as f64])
        .collect::<Vec<_>>();
    let curve_plot = poloto::build::line("curve", curve_points);

    let derivative_1_points = range
        .clone()
        .map(|x| [x as f64, curve.derivative(1, &x, &mut buffer).real as f64])
        .collect::<Vec<_>>();
    let derivative_1_plot = poloto::build::line("derivative 1", derivative_1_points);

    let derivative_2_points = range
        .clone()
        .map(|x| [x as f64, curve.derivative(2, &x, &mut buffer).real as f64])
        .collect::<Vec<_>>();
    let derivative_2_plot = poloto::build::line("derivative 2", derivative_2_points);

    let derivative_3_points = range
        .clone()
        .map(|x| [x as f64, curve.derivative(3, &x, &mut buffer).real as f64])
        .collect::<Vec<_>>();
    let derivative_3_plot = poloto::build::line("derivative 3", derivative_3_points);

    let derivative_4_points = range
        .clone()
        .map(|x| [x as f64, curve.derivative(4, &x, &mut buffer).real as f64])
        .collect::<Vec<_>>();
    let derivative_4_plot = poloto::build::line("derivative 4", derivative_4_points);

    let plots = plots!(
        knots_plot,
        curve_plot,
        derivative_1_plot,
        derivative_2_plot,
        derivative_3_plot,
        derivative_4_plot
    );

    let plotter = simple_fmt!(plots, title_name, "x", "y");

    // Write the plot to a file
    let folder = Path::new("plots");
    let file = folder.join(file_name);

    std::fs::create_dir_all(folder).expect("Failed to create the folder for the plots");

    let mut file = File::create(file).expect("Failed to create the file for plot");

    write!(file, "{}", poloto::disp(|f| plotter.simple_theme(f)))
        .expect("Failed to create write the plot");
}

fn create_bspline_example(
    degree: usize,
    knots: Vec<f32>,
    control_points: Vec<f32>,
) -> BSplineCurve<f32, NoInfinitesimal> {
    let dim = Dim(0);
    let jet_control_points = control_points
        .into_iter()
        .map(|p| Jet::new(p, NoInfinitesimal))
        .collect();
    let bspline = BSpline::<f32>::new(degree, knots).unwrap();
    BSplineCurve::new(dim, bspline, jet_control_points).unwrap()
}

fn bspline_example_0_a() -> BSplineCurve<f32, NoInfinitesimal> {
    create_bspline_example(0, vec![1.0], vec![2.0])
}

fn bspline_example_0_b() -> BSplineCurve<f32, NoInfinitesimal> {
    create_bspline_example(0, vec![1.0, 2.0, 3.0], vec![2.0, 3.0, 1.0])
}

fn bspline_example_1_a() -> BSplineCurve<f32, NoInfinitesimal> {
    create_bspline_example(1, vec![1.0, 2.0], vec![1.0, 3.0])
}

fn bspline_example_1_b() -> BSplineCurve<f32, NoInfinitesimal> {
    create_bspline_example(1, vec![1.0, 2.0, 3.0, 4.0], vec![1.0, 3.0, 2.0, 5.0])
}

fn bspline_example_2_a() -> BSplineCurve<f32, NoInfinitesimal> {
    create_bspline_example(2, vec![1.0, 2.0, 3.0], vec![1.0, 3.0, 2.0, 4.0])
}

fn bspline_example_2_b() -> BSplineCurve<f32, NoInfinitesimal> {
    create_bspline_example(
        2,
        vec![1.0, 2.0, 3.0, 4.0, 5.0],
        vec![1.0, 3.0, 2.0, 4.0, 1.0, 5.0],
    )
}

fn bspline_example_3_a() -> BSplineCurve<f32, NoInfinitesimal> {
    create_bspline_example(
        3,
        vec![1.0, 2.0, 3.0, 4.0],
        vec![1.0, 3.0, 2.0, 4.0, 1.0, 5.0],
    )
}

fn bspline_example_3_b() -> BSplineCurve<f32, NoInfinitesimal> {
    create_bspline_example(
        3,
        vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        vec![1.0, 3.0, 2.0, 4.0, 1.0, 5.0, 2.0, 3.0],
    )
}
