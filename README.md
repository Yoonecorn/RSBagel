# RSBagel

RSBagel is a Rust implementation of BAGEL2 for Bayesian analysis of pooled CRISPR fitness screens.

**Supported source platform: Apple Silicon macOS (`aarch64-apple-darwin`).**

## Documentation

- [Quick Start](docs/RSBagel-Quick-Start.md): begin with the synthetic example.
- [User Manual](docs/RSBagel-User-Manual.md): routine analysis and interpretation.
- [Advanced CLI Reference](docs/RSBagel-Advanced-CLI-Reference.md): command options and profiles.
- [Developer Integration Guide](docs/RSBagel-Developer-Integration-Guide.md): Rust and C integration.
- [Manuscript](docs/RSBagel-Manuscript.pdf): methods, comparisons and limitations.
- [Beginner example](examples/beginner/README.md) and [third-party notices](docs/BAGEL2-LICENSE.txt).

## What RSBagel does

RSBagel calculates fold changes, gene Bayes factors and precision–recall results from guide count tables and reference gene sets. It provides a command-line interface, a Rust core and a C interface. Individual sgRNAs in a Brunello-style screen are separate observations; a four-guide single-target Inzolia-style array is one observation and must not be split into four artificial observations.

The default NumPy-compatible sampling profile follows the pinned BAGEL2 compatibility path. Native ChaCha8 sampling is deterministic and experimental; its stream does not reproduce NumPy draws. Observable agreement does not imply universal internal bitwise identity.

## Requirements and installation

Use an Apple Silicon Mac with Rust and Cargo installed. From a checkout, build the executable:

```sh
unset CARGO_TARGET_DIR
cargo build --release --locked --target aarch64-apple-darwin
```

## Quick Start

Follow the [Quick Start](docs/RSBagel-Quick-Start.md) for setup and the complete synthetic example. To inspect the executable from the repository root:

```sh
target/aarch64-apple-darwin/release/rsbagel --help
```

The example produces `results.foldchange`, `results.normed_readcount`, `results.bf`, and `results.pr`, plus metadata sidecars. The manual describes the columns and their interpretation.

## Current limits

Qualification applies to the tested macOS ARM64 source configuration and described datasets. Other platforms and the native sampling profile are not qualified. Scientific results depend on input quality, control sets and study design; computational agreement does not establish biological validity. Memory use and performance depend on workload.

## Citation and license

Cite RSBagel by repository and exact commit, and cite BAGEL2 as described in the [manuscript](docs/RSBagel-Manuscript.pdf). No DOI is assigned here. RSBagel is licensed under the [MIT License](LICENSE); see the [BAGEL2](docs/BAGEL2-LICENSE.txt) and [NumPy](docs/NUMPY-LICENSE.txt) notices.
