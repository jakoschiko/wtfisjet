# wtfisjet

Provides an implementation for jet and some jet-based utilities.

## Status of this crate

The author does not consider this crate as stable yet. However, it is well documented and
tested, so give it a try!

## Nightly rust

This crate uses some nightly features. Hence it's necessary to use the nightly rust compiler.
Hopefully this won't be necessary anymore in the future.

Nightly features used by this crate:
- `generic_associated_types` is used for the inner [`Iterator`] types of [`Infinitesimal`].
We can expect the stabilization in the near future, see [The push for GATs stabilization].
- `total_cmp` allows us to use function like [`slice::binary_search_by`] for [`f32`] and
[`f64`]. Although it's handy, this could also be solved with workarounds. See the
[tracking issue for `total_cmp`].

[The push for GATs stabilization]: https://blog.rust-lang.org/2021/08/03/GATs-stabilization-push.html
[tracking issue for `total_cmp`]: https://github.com/rust-lang/rust/issues/72599

## Features

These features are available:

- Implementation for jets that can be used for AD.
- Multiple implementation for the infinitesimal part of the jet with different performance
characteristics.
- Implementation for B-splines that uses jets as control points.
- Implementation for Newton's method that uses jets for calculating the derivative.
- Many feature flags that allows to include only the features you need.

## WTF is jet?

A jet is a magic number. If you use jets for your calculation instead of ordinary numbers,
you get not only the result but also the derivatives with respect to the inputs. This is
called [automatic differentiation (AD)]. There are different ways of accomplishing AD,
[dual numbers] are one of them and a jet is an n-dimensional dual number. The [ceres library]
has a good explanation of how all of this works.

[automatic differentiation (AD)]: https://en.wikipedia.org/wiki/Automatic_differentiation
[dual numbers]: https://en.wikipedia.org/wiki/Dual_number
[ceres library]: http://ceres-solver.org/automatic_derivatives.html#dual-numbers-jets

This library implements jet with the struct [`Jet`]. The most important thing you need to know
about [`Jet`] is that it consists of two parts:
- The real part is the actual result.
- The infinitesimal part is the magic part that calculates the derivatives. It's represented
by the trait [`Infinitesimal`]. This crate provides different implementations with different
performance characteristics.

So, how can you use [`Jet`]? Instead of writing this:

```rust
fn foo(x: f32) -> f32 {
    x * x - 1.0
}

let result = foo(2.0);
println!("{}", result);
// Output: 3
```

You write this:

```rust
use wtfisjet::{DenseInfinitesimal, Dim, Infinitesimal, Jet};

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
```

As mentioned before, jets are n-dimensional. That means you can calculate the derivatives with
respect to multiple inputs:

```rust
use wtfisjet::{DenseInfinitesimal, Dim, Infinitesimal, Jet};

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
```

Please note that functions `foo` and `bar` in the previous examples don't hide [`Jet`]
behind a type parameter. Of course you could do this:

```rust
use std::ops::Mul;
use wtfisjet::{DenseInfinitesimal, Dim, Infinitesimal, Jet};

// Note: Not recommended!
fn baz<N: Clone + Mul<Output=N>>(x: N) -> N {
    x.clone() * x
}

type MyJet = Jet<f32, DenseInfinitesimal<f32>>;
let dim = Dim(1);

let result = baz(MyJet::variable(2.0, 0, dim));
```

However, this is not recommended. [`Jet`] has many pitfalls and it's helpful to see
that your calculation is dealing with jets and not ordinary numbers. Additionally you can
often improve the performance by using special functions that reduce the number of
necessary [`Infinitesimal`] operations (e.g. you can use `x.square()` instead of
`x.clone() * x`).

If you're not interested in the derivatives but want to use a [`Jet`] based function,
you can use [`NoInfinitesimal`]:

```rust
use wtfisjet::{Dim, Infinitesimal, Jet, NoInfinitesimal};

fn foo<I: Infinitesimal<f32>>(x: Jet<f32, I>) -> Jet<f32, I> {
    x.clone() * x - 1.0
}

type MyJet = Jet<f32, NoInfinitesimal>; // NoInfinitesimal is a zero sized type (ZST)
let dim = Dim(0); // NoInfinitesimal requires zero dimensions

let result = foo(MyJet::constant(2.0, dim)).real;
println!("{}", result);
// Output: 3
```

Assuming that Rust's zero-overhead principle is not a lie, this will be as fast as using
a `f32` based function.

## Example: Curve fitting

What is [curve fitting]? Imagine you have a tool that takes some variables and produce a curve.
You want to produce a curve that adheres some constraints. How do you choose the variables?
Curve fitting is the process of finding those variables.

This library provides utilities for simple curve fitting:
- You can use B-splines as a tool for constructing curves. The variables are the control points
of the B-spline. The constraints are points that should be matched by the B-spline curve.
- You can use Newton's method to approximate a solution for the variables.
- You can use jets to calculate the derivative that is necessary for Newton's method.

The example [curve_fitting.rs] shows how this works. It produces the following curve:

![curve fitting result]

[curve fitting]: https://en.wikipedia.org/wiki/Curve_fitting
[curve_fitting.rs]: https://github.com/jakoschiko/wtfisjet/blob/master/examples/curve_fitting.rs
[curve fitting result]: https://raw.githubusercontent.com/jakoschiko/wtfisjet/master/plots/curve_fitting.svg

## Feature flags

These feature flags can be used to customize the crate.

### `big-rational-number` (enabled by default)

If enabled, the type `num-rational::BigRational` will implement [`Number`].
This is very useful for writing tests with arbitrary precision. It will
require the dependency `num-rational`.

### `sparse-infinitesimal` (enabled by default)

If enabled, the type `wtfisjet::SparseInfinitesimal` is available. It will
require the dependency `intmap`.

### `const-infinitesimal` (enabled by default)

If enabled, the type `wtfisjet::ConstInfinitesimal` is available. It will
require the dependency `array-init`.

### `b-spline` (enabled by default)

If enabled, the module `wtfisjet::bspline` is available.

### `newton` (enabled by default)

If enabled, the module `wtfisjet::newton` is available. It will require the dependency
`rulinalg`.

### `dice` (disabled by default)

If enabled, the module `wtfisjet::dice` is available. It provides generators for
random test data. It will require the dependency `dicetest`.

### `asserts` (disabled by default)

If enabled, the module `wtfisjet::asserts` is available. It provides assertions for tests.
It will require the dependency `dicetest` and `diceprop`.

## License

Licensed under either of

 * Apache License, Version 2.0
   ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
 * MIT license
   ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

at your option.

## Contribution

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in the work by you, as defined in the Apache-2.0 license, shall be
dual licensed as above, without any additional terms or conditions.
