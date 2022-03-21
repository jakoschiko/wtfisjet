#!/bin/bash
set -e

echo "STEP wtfisjet: cargo +nightly fmt"
cargo +nightly fmt

echo "STEP wtfisjet: cargo +nightly build"
cargo +nightly build

echo "STEP wtfisjet: cargo +nightly build --no-default-features"
cargo +nightly build --no-default-features

echo "STEP wtfisjet: cargo +nightly build --all-features"
cargo +nightly build --all-features

echo "STEP wtfisjet: cargo +nightly test -- --format=terse"
cargo +nightly test -- --format=terse

echo "STEP wtfisjet: cargo run --example bspline_examples"
cargo run --example bspline_examples

echo "STEP wtfisjet: cargo +nightly clippy"
cargo +nightly clippy

echo "STEP wtfisjet: cargo +nightly doc --no-deps"
cargo +nightly doc --no-deps

echo "STEP wtfisjet: cargo readme > README.md"
cargo readme > README.md
