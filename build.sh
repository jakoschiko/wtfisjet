#!/bin/bash
set -e

echo "STEP wtfisjet: cargo fmt"
cargo fmt

echo "STEP wtfisjet: cargo build"
cargo build

echo "STEP wtfisjet: cargo build --no-default-features"
cargo build --no-default-features

echo "STEP wtfisjet: cargo build --all-features"
cargo build --all-features

echo "STEP wtfisjet: cargo test -- --format=terse"
cargo test -- --format=terse

echo "STEP wtfisjet: cargo run --example bspline_examples"
cargo run --example bspline_examples

echo "STEP wtfisjet: cargo run --example curve_fitting"
cargo run --example curve_fitting

echo "STEP wtfisjet: cargo clippy"
cargo clippy

echo "STEP wtfisjet: cargo doc --no-deps"
cargo doc --no-deps

echo "STEP wtfisjet: cargo readme > README.md"
cargo readme > README.md
