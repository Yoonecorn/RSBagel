# Contributing

Use an Apple Silicon Mac with Git, Rust/Cargo and Apple command-line developer tools. Clone the repository and build as described in the [README](README.md). From the repository root run:

```sh
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --locked -- -D warnings
cargo test --workspace --locked
cargo build --release --locked --target aarch64-apple-darwin
python3 scripts/check_doc_links.py
```

Follow the [Quick Start](docs/RSBagel-Quick-Start.md) to build and run the synthetic beginner example and compare all four tables. Describe any scientific-output change with reproducible inputs, commands, expected and observed results, and its effect on sampling, precision, ordering or interpretation. Do not upload sensitive or restricted data.

Changes to the scientific contract or outputs require explicit qualification. Publication benchmarks, external validation, and scientific qualification gates are not routine pull-request tests. See the [developer guide](developer/RSBagel-Developer-Integration-Guide.md) for ignored external-data and qualification tests.
