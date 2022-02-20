# wtfisjet

Provides an implementation for jets with some utilities.

A jet is n-dimensional dual number that can be used for automatic derivation.
See the [ceres library] for more information.

[ceres library]: http://ceres-solver.org/automatic_derivatives.html#dual-numbers-jets

## Feature flags

These feature flags can be used to customize the crate.

### `big-rational` (disabled by default)

If enabled, the type `num-rational::BigRational` will implement [`Number`].
This is very useful for writing tests with arbitrary precision. It will
require the dependency `num-rational`.

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
