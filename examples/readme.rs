// Contains the examples used in the readme.

use wtfisjet::{DenseInfinitesimal, Dim, Infinitesimal, Jet, NoInfinitesimal};

fn main() {
    {
        println!("# foo with f32");

        fn foo(x: f32) -> f32 {
            x * x - 1.0
        }

        let result = foo(2.0);
        println!("{}", result);
        // Output: 3

        println!()
    }

    {
        println!("# foo with Jet");

        fn foo<I: Infinitesimal<f32>>(x: Jet<f32, I>) -> Jet<f32, I> {
            x.clone() * x - 1.0
        }

        type MyJet = Jet<f32, DenseInfinitesimal<f32>>;
        let dim = Dim(1);

        let result = foo(MyJet::variable(2.0, 0, dim));
        println!("{}", result);
        // Output: (3 + [4]h)
        //          ^    ^
        //          |    |
        //          |   derivative of foo with respect to x
        //          |
        //         result of foo

        println!()
    }

    {
        println!("# bar with Jet and two variables");

        fn bar<I: Infinitesimal<f32>>(x: Jet<f32, I>, y: Jet<f32, I>) -> Jet<f32, I> {
            (x.clone() * x - 1.0) / y
        }

        type MyJet = Jet<f32, DenseInfinitesimal<f32>>;
        let dim = Dim(2); // Because we have now two variables, we need at least 2 dimensions

        let result = bar(
            // The different variables uses different indices
            MyJet::variable(2.0, 0, dim),
            MyJet::variable(4.0, 1, dim),
        );
        println!("{}", result);
        // Output: (0.75 + [1, -0.1875]h)
        //          ^^^^    ^  ^^^^^^^
        //           |      |     |
        //           |      |    derivative of bar with respect to y
        //           |      |
        //           |     derivative of bar with respect to x
        //           |
        //         result of bar

        println!()
    }

    {
        println!("# foo with Jet and NoInfinitesimal");

        fn foo<I: Infinitesimal<f32>>(x: Jet<f32, I>) -> Jet<f32, I> {
            x.clone() * x - 1.0
        }

        type MyJet = Jet<f32, NoInfinitesimal>; // NoInfinitesimal is a zero sized type (ZST)
        let dim = Dim(0); // NoInfinitesimal requires zero dimensions

        let result = foo(MyJet::constant(2.0, dim)).real;
        println!("{}", result);
        // Output: 3

        println!()
    }
}
