# RSBagel example

This is an offline teaching package for macOS on Apple Silicon (arm64).
Read the [Quick Start](../../docs/RSBagel-Quick-Start.md) before running it.
The executable is a locally built RSBagel 0.1.0 program, not a public release.
Inputs are synthetic, not patient or experimental data.

The three tab-separated input files contain 360 independent reagent observations
for 120 artificial genes, with baseline T0 and endpoint samples R1/R2.
G000–G039 are simulated essential controls; G040–G079 are simulated nonessential
controls. G080–G119 are simulated intermediate-depletion genes.
Each gene has three separate reagent observations. These are individual-guide
style inputs, not four-guide arrays. R1/R2 are synthetic endpoint samples;
the example does not demonstrate real biological replication or validate hits.

Inputs are exact copies from tests/golden/synthetic-bootstrap-final/inputs in
the source repository at commit 117abd6931f3e0528f6afd9aa6b7a60e5d7ef8db.
The inputs are synthetic and are included for learning the workflow.

The run script uses NumPy-compatible sampling, seed 42, one worker and default
10-fold cross-validation. No Python or NumPy installation is needed.
It refuses to overwrite an existing outputs folder. Its comparison checks all
four scientific tables byte for byte against expected/. Metadata sidecars
record configuration, command arguments, paths, input hashes, and build identity;
their bytes can differ when folders or commands change. The metadata sidecars are not byte compared.

Build and copy the executable into bin/ as described in the Quick Start. licenses/ contains retained notices. This local packaging does not
extend platform or scientific qualification.
