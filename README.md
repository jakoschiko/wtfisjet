# wtfisjet

Provides an implementation for jets with some utilities.

A jet is n-dimensional dual number that can be used for automatic derivation.
See the [ceres library] for more information.

[ceres library]: http://ceres-solver.org/automatic_derivatives.html#dual-numbers-jets

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
