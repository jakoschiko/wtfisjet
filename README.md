# wtfisjet

Provides an implementation for jets with some utilities.

A jet is n-dimensional dual number that can be used for automatic derivation.
See the [ceres library] for more information.

[ceres library]: http://ceres-solver.org/automatic_derivatives.html#dual-numbers-jets

## Nightly rust

This crate uses the nightly feature `generic_associated_types`. Hence it's necessary to use
the nightly rust compiler. Hopefully this won't be necessary anymore in the near future,
see [The push for GATs stabilization].

[The push for GATs stabilization]: https://blog.rust-lang.org/2021/08/03/GATs-stabilization-push.html

## Feature flags

These feature flags can be used to customize the crate.

### `big-rational-number` (disabled by default)

If enabled, the type `num-rational::BigRational` will implement [`Number`].
This is very useful for writing tests with arbitrary precision. It will
require the dependency `num-rational`.

### `sparse-infinitesimal` (disabled by default)

If enabled, the type `wtfisjet::SparseInfinitesimal` is available. It will
require the dependency `intmap`.

### `const-infinitesimal` (disabled by default)

If enabled, the type `wtfisjet::ConstInfinitesimal` is available. It will
require the dependency `array-init`.

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
