# RSBagel Developer Integration Guide

RSBagel 0.1.0. **Wet-lab users can skip this guide.** Use the
[User Manual](../docs/RSBagel-User-Manual.md) for the routine workflow.
Use the executable’s `--help` option for command-line options.
A repository is the source-code folder. A binary is a compiled program or library.

These interfaces are for software integrators. Wet-lab users can skip this section.

### Rust library

`rsbagel-core` is the in-memory Rust library. Its main entry points are `run_fc`, `run_bf`, and `run_pr`. `Matrix` contains reagent IDs, gene labels, column labels, and a flat row-major `f64` value buffer. Configuration fields are documented in the source types. Library numeric-column indices begin at **zero**, unlike the CLI's routine one-based selectors.

`FcConfig::default()` selects numeric column zero, minimum reads zero, and pseudocount five. `BfConfig::default()` selects NumPy-compatible sampling, seed 42, ten-fold cross-validation, 1,000 bootstrap iterations when that mode is selected, the automatic worker policy, and the legacy network profile. The CLI's omitted seed instead uses time; do not assume the two entry points have the same seed default.

For CLI-compatible composition, use `fc_to_bf_input` between FC and BF, then `bf_to_pr_input` between BF and PR. These bridges preserve the source three-decimal boundaries and the alphabetical BF serialization boundary without writing files. Passing raw full-precision stage results directly is a different numerical path. It must not be described as automatically identical to the text-file CLI chain.

### C ABI

A C ABI, or application binary interface, describes how compiled software calls functions in a library. RSBagel's ABI version is 1. The header is `crates/rsbagel-ffi/include/rsbagel.h` in the source repository. The library offers FC, BF, and PR entry points, result accessors, progress callbacks, and cancellation tokens. It performs no disk I/O or subprocess analysis.

Strings are UTF-8 pointer-and-length views, not NUL-terminated strings. Matrices use row-major `double` values. Inputs are borrowed synchronously and copied into the owned core matrix. Input allocations must remain valid for the call. Result handles own their returned buffers. Check status before reading results. Free each returned handle exactly once, including failures. Borrowed result views expire when their handle is freed.

| Status | Meaning |
|---|---|
| 0 | Success |
| 1 | Invalid argument or input |
| 2 | Numerical model error |
| 3 | Preserved legacy failure |
| 4 | Cancelled |
| 5 | I/O error |
| 6 | Internal error |

Empty configuration bytes select defaults. Numeric column indices in JSON begin at zero. Reference inputs are arrays of actual gene strings without file headers. The host must retain its complete configuration and input identities.

FC buffers zero and one are raw fold changes and normalized counts. Optional buffer two is `cli_bf_input`. BF buffers zero to three are gene summaries, guide means, guide standard deviations, and guide totals. Optional BF buffer four is `cli_pr_input`; use its accompanying `cli_pr_genes` labels. Check availability metadata before using either optional ready-input buffer. Do not apply replicate scaling a second time to a ready buffer.

Raw BF gene fields are `raw_mean`, `mean`, `std`, `normalized`, `boosted`, `boosted_std`, and `boosted_normalized`. Raw buffers preserve IEEE floating-point nonfinite values. JSON may represent them as `null`. Metadata distinguishes absent optional fields from present nonfinite values. Serialization scaling is field-specific. PR returns score, recall, precision, and FDR in core ranked order.

Progress callbacks can run on worker threads. Return nonzero to cancel. Cancellation is cooperative. A cancellation token must outlive the active call and any concurrent request. Do not free a result while another thread reads it. These interfaces do not establish a complete graphical application or a SwiftUI product. Existing C and Swift interface checks do not expand scientific or platform qualification.

## Matrix shapes and buffer layouts

An input matrix has one gene label per reagent ID. Its value buffer length is
row count multiplied by numeric column count, in row-major order. Rows follow
reagent IDs and columns follow numeric column labels. C views use `double` for
Rust `f64` and `size_t` for lengths. Host architecture and pointer width must match.

| Workflow | Buffer | Shape and meaning |
|---|---|---|
| FC | 0 `fold_change` | Retained rows × noncontrol sample columns |
| FC | 1 `normalized` | Retained rows × sample/combined-control columns |
| FC | Optional 2 `cli_bf_input` | FC rounded at the CLI three-decimal boundary for BF |
| BF | 0 `genes` | First-seen genes × seven fields, in the order listed above |
| BF | 1 `guide_means` | Reagents × selected replicate columns |
| BF | 2 `guide_stds` | Reagents × selected replicate columns |
| BF | 3 `guide_totals` | Reagents × one total score column |
| BF | Optional 4 `cli_pr_input` | Ready scores paired with metadata `cli_pr_genes` |
| PR | 0 `precision_recall` | Ranked genes × score, recall, precision, FDR |

Use metadata for labels, counts, optional-field presence, and per-field
serialization scaling. Raw result ordering is not automatically CLI alphabetical
ordering. Keep raw buffers separate from ready-input buffers. An unavailable
bridge produces availability/error metadata; it must not be treated as a valid
empty successful input.

## Pointer ownership, thread safety, and lifetimes

Nonempty pointers must be aligned, readable, live allocations of the declared
lengths. Null pointers are permitted only for zero-length input views. Validation
of nullness, alignment, UTF-8 and dimensions cannot prove that an arbitrary
nonnull C pointer points to valid memory. That remains the host's responsibility.
Do not separately free or mutate borrowed result views. Check status before
interpreting buffer counts or numeric values. Read-only accessors can be used
concurrently while their owned handle stays live. A handle must not be freed
while another thread reads it.

A foreign callback must not unwind across the C boundary. Its context must be
safe for worker-thread access. It must not free in-use input buffers or tokens.
The cancellation token must outlive the workflow and all concurrent requests.
Reset a token only before a later run, not during an active call. Free it after
all calls and requests have completed. Cancellation occurs at cooperative progress
boundaries; it is not forced thread termination. Contained workflow panics do
not make arbitrary memory corruption or process-aborting failure recoverable.

## Qualification and technical environment

The documented evidence platform is macOS arm64 / Rust 1.98 / Python 3.11.
These interface contracts do not establish another platform or host UI.
The supplied executable needs no Rust, Python, or NumPy installation.
NumPy-compatible sampling is the default for CLI, core and omitted C JSON fields.
Explicit native ChaCha8 is implemented and deterministic but remains experimental,
with the qualification boundary explained in the advanced reference.

## External HAP1 regression tests

The two HAP1 regression tests are retained but ignored in an ordinary checkout because the approximately 56 MB external HAP1 golden dataset is intentionally excluded. They are not part of the clean-checkout passing test count. An authorized developer with the unchanged golden dataset can place it at `tests/golden/hap1-cv-final/` and the matching pinned BAGEL2 reference gene lists at `upstream/BAGEL2/CEGv2.txt` and `upstream/BAGEL2/NEGv1.txt`, then run from the repository root:

```sh
cargo test --locked -p rsbagel-core bf::tests::hap1_replay_stages -- --ignored --exact
cargo test --locked -p rsbagel-core resampling::tests::hap1_numpy_plan_matches_reference -- --ignored --exact
```

These commands require the external fixtures; they do not download or regenerate them. The 1,000-iteration native bootstrap distributional gate is a bounded scientific qualification gate, not a routine software test. It is ignored in ordinary testing and was not run during this assembly. Only during an explicitly approved future qualification campaign, an authorized developer with the unchanged external qualification fixture at `tests/external/qualification/bootstrap1000/` can run:

```sh
cargo test --locked -p rsbagel-core bf::tests::native_bootstrap_distributional_gate -- --ignored --exact
```

The qualification fixture is not distributed in this snapshot. Four additional evidence-generation and diagnostic tests were already ignored in the original suite; they are not run by ordinary testing.
