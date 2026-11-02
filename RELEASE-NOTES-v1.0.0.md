# RSBagel v1.0.0 release candidate

RSBagel implements the BAGEL2 pooled CRISPR fitness-screen workflow. Source use is qualified for the evaluated Apple Silicon macOS configuration; other operating systems are not qualified. There is no binary distribution. Build from source using the [README](README.md).

The `fc`, `bf`, and `pr` stages calculate fold changes, estimate gene Bayes factors, and evaluate ranked genes against essential and nonessential reference sets. The engine provides a CLI, Rust API, and C ABI. NumPy-compatible sampling is the default compatibility profile. Native ChaCha8 is deterministic and experimental; replay consumes a supplied partition plan. Equal numerical seeds do not imply equal streams across generators.

Run the synthetic [beginner example](examples/beginner/README.md) via the [Quick Start](docs/RSBagel-Quick-Start.md). See the [User Manual](docs/RSBagel-User-Manual.md), [developer guide](developer/RSBagel-Developer-Integration-Guide.md), and [manuscript](docs/RSBagel-Manuscript.pdf).

Individual sgRNAs are separate observations; a four-guide single-target array is one observation. Computational agreement does not prove biological validity. Results and performance are bounded by the evaluated data and configurations; universal bitwise identity, speedup or memory reduction is not claimed.

Routine tests do not establish scientific qualification. Two HAP1 regression tests are ignored because their external golden dataset is not distributed; the native 1,000-iteration bootstrap test is an ignored scientific qualification gate. The four additional evidence-generation or diagnostic tests are also ignored. Ignored tests are not counted as passing.

RSBagel uses the [MIT license](LICENSE); [BAGEL2 and NumPy notices](third-party-licenses/) are included.
