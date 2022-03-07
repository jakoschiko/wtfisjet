use std::fmt::Debug;

use float_eq::float_eq;

use crate::{Infinitesimal, Jet};

pub fn assert_jet_eq<I: Infinitesimal<f64> + Debug>(
    left: Jet<f64, I>,
    right: Jet<f64, I>,
    tol: f64,
) {
    assert_eq!(
        left.infinitesimal.dim(),
        right.infinitesimal.dim(),
        "Dimension counts the jets are not equal\n left: {:?}\nright: {:?}",
        left.infinitesimal.dim(),
        right.infinitesimal.dim(),
    );

    assert!(
        float_eq!(left.real, right.real, abs <= tol),
        "Real parts of the jets are not equal\n  tol: {}\n diff: {}\n left: {:?}\nright: {:?}",
        tol,
        (left.real - right.real).abs(),
        left,
        right,
    );

    for (idx, (&left_i, &right_i)) in left
        .infinitesimal
        .dense_elems()
        .zip(right.infinitesimal.dense_elems())
        .enumerate()
    {
        assert!(
            float_eq!(left_i, right_i, abs <= tol),
            "Infinitesimal parts of the jets are not equal at dimension {}\n  tol: {}\n diff: {}\n left: {:?}\nright: {:?}",
            idx,
            tol,
            (left_i - right_i).abs(),
            left,
            right,
        )
    }
}
