# RSBagel Advanced CLI Reference

RSBagel 0.1.0. CLI means command-line interface: the text interface used in Terminal.
Read the [User Manual](RSBagel-User-Manual.md) for routine work.
Use the [Developer Integration Guide](RSBagel-Developer-Integration-Guide.md) for library calls.

This guide preserves the previously validated option details and command text.
Advanced legacy combinations below are **not recommended for routine use**.
A parser accepting a setting does not establish its scientific validity.

## Technical platform note

The documented validation environment is macOS arm64 / Rust 1.98 / Python 3.11.
The local executable check used macOS 27.0.1 arm64 and Rust/cargo 1.98.0.
Python describes the reference environment; the supplied executable does not
require Python or NumPy. Other platforms are not qualified by that evidence.
This local package is not a signed or notarized public release.

## Complete options, aliases, exit statuses, and failures

The option tables below were checked against the executable and small synthetic runs. In these tables, “required” means the option must be supplied. A file option takes a file name. The spelling includes two leading hyphens for long names and one for short names. Optional settings should be changed only when you understand their effect.

### Top-level commands

`fc` calculates log2 fold changes from read counts. `bf` calculates Bayes-factor scores. `pr` ranks gene scores and compares reference labels. `version` prints BAGEL2 compatibility version 2.0/build 115 and RSBagel 0.1.0. `--version` or `-V` prints RSBagel 0.1.0 alone. `--help` or `-h` shows help. The generated `help` command also displays subcommand help. A missing required option, unknown option, or invalid option syntax exits with status 2. A data or computation failure exits with status 1. Success exits with status 0.

### fc

| Option | Default / accepted input | Effect / failure |
|---|---|---|
| `-i`, `--read-count-file` | Required existing tab-separated file | Count matrix with reagent ID, gene ID, then numeric columns. Missing file fails. |
| `-o`, `--output-label` | Required output prefix | Writes prefix `.foldchange`, `.normed_readcount`, `.metadata.json`. Missing parent folder fails. |
| `-c`, `--control-columns` | Required comma-separated names or numeric positions | First numeric column is 1. Control columns are excluded from `.foldchange`. Names must match exactly. Zero or nonexistent controls fail. Negative indexing follows legacy semantics; use names or positive positions for clarity. |
| `-m`, `--min-reads` | Signed integer; default 0 | Filters rows using summed raw control reads plus one pseudocount, before normalization (not the mean). |
| `--pseudo-count` | Signed integer; default 5 | Pseudocount used for FC normalization/calculation. Use positive counts for a normal analysis. Legacy alias `-Np`. |
| `--dump-json` | Optional output JSON path | Writes raw FC input/config/result and provenance. Missing parent folder fails. |
| `-h`, `--help` | Flag | Prints help, no analysis. |

### bf

| Option | Default / accepted input | Effect / failure |
|---|---|---|
| `-i`, `--fold-change` | Required existing TSV | Usually `.foldchange` from `fc`. Duplicate reagent IDs fail. |
| `-o`, `--output-file` | Required file path | Writes exactly this file plus `.metadata.json`; no `.bf` suffix is added automatically. |
| `-e`, `--essential-genes` | Required existing file | One header line, then gene names in first tab-separated field. First line is always skipped. |
| `-n`, `--non-essential-genes` | Required existing file | Same format as essentials. Training must contain usable, varied reference observations. |
| `-c`, `--columns-to-test` | Required comma-separated names or positions | First numeric FC column is 1. Zero, out-of-range, or no matching names fails. Negative `-1` fails; avoid negative positions. |
| `-w`, `--network-file` | Optional existing file | Headerless tab-separated gene-pair network; adds network scores for gene output. Guide output skips network correction. |
| `-m`, `--filter-multi-target` | Flag; off | Requires `--align-info`. Filters/corrects according to alignment data. |
| `--loci-without-mismatch` | Nonnegative integer; 10 | Maximum perfect-match locus count. Used with multi-target filtering. Legacy alias `-m0`. |
| `--loci-with-mismatch` | Nonnegative integer; 10 | Maximum mismatch-locus count. Used with filtering. Legacy alias `-m1`. |
| `--align-info` | Optional existing file | Headerless five-field alignment information. Needed with `-m`. Even without `-m`, a supplied nonexistent path fails while recording metadata. |
| `-b`, `--use-bootstrapping` | Flag; off | Bootstrap resampling instead of default cross-validation. |
| `--no-resampling` | Flag; off | Uses one fit with every gene in both training and scoring sets. Overrides iteration generation. Legacy alias `-NS`. This is not independent validation. |
| `--use-small-sample` | Flag; off | 100 reference-gene bootstrap fits; overrides bootstrap/CV choice. Unlabelled genes are scored each iteration. |
| `-N`, `--no-of-cross-validations` | Nonnegative integer; 10 | Default cross-validation fold count. Zero fails during CV. Too many folds for the gene universe fails. |
| `--bootstrap-iterations` | Nonnegative integer; 1000 | Iteration count with bootstrap. Zero fails when bootstrap iteration generation is active. Legacy alias `-NB`. |
| `-r`, `--sgrna-bayes-factors` | Flag; off | Produces reagent-level scores instead of gene-level rows. Do not feed these rows directly into gene-level PR. |
| `-f`, `--equalise-sgrna-no` | Optional signed integer; unset | Adds guide-count-normalized gene score `NormBF`. Zero disables normalization. Negative integers parse and run but have no recommended scientific use. |
| `-p`, `--equalise-rep-no` | Optional signed integer; unset | Scales designated BF output scores to requested replicate count. Zero disables. Negative values parse and run but have no recommended scientific use. |
| `-s`, `--seed` | Unsigned integer; clock-derived if absent | NumPy mode accepts 0–4294967295. Native accepts 0–18446744073709551615. Negative and noninteger values fail parsing. Above NumPy range fails during computation. Record an explicit seed. |
| `-t`, `--run-test-mode` | Flag; off | Adds `.traininfo`. With network gene output also writes `.netscore`. This is a diagnostic flag, not an automated validation claim. |
| `--rng-compat` | `numpy` (default), `rust` | `numpy` is production default. `rust` selects experimental, unqualified native ChaCha8. Other strings fail. |
| `--network-profile` | `legacy` (default), `corrected` | Controls optional network semantics. Corrected separates boosted and normalized fields. Other strings fail. |
| `--replay-manifest` | Optional existing JSON | Reuses explicit partitions. Gene identities/order, contiguous iteration IDs, complete disjoint membership, and requested iteration count must match. Incorrect JSON/membership fails. |
| `--workers` | `auto` (default) or positive integer | Resolved automatic count is recorded in metadata. Zero, negative, fractional or text other than `auto` fails. Explicit counts have no automatic cap. |
| `--threads` | Alias for `--workers` | Same behavior. |
| `-v`, `--verbose` | Flag; off | Prints resolved worker count to standard error. |
| `--compat` | Flag; off | Explicit alias for NumPy mode. Takes precedence even when `--rng-compat rust` is also present. |
| `--dump-json` | Optional JSON output path | Writes raw result/config/stage details and provenance. Can be much larger than TSV. |
| `-h`, `--help` | Flag | Prints help. |

Repeated value options use the last value. The parser does not reject combined resampling flags. Selection priority is small-sample, then bootstrap, then CV. `--no-resampling` overrides partition generation but does not reset the selected output schema. For clear analysis records choose one resampling method. The CLI metadata `ranking_column` recommends a network score; `pr` still defaults to `BF`, so select the desired network column explicitly.

### pr

| Option | Default / accepted input | Effect / failure |
|---|---|---|
| `-i`, `--bayes-factors` | Required existing TSV | Gene-level BF table with header and numeric score column. |
| `-o`, `--output-file` | Required output path | Writes exactly this file and `.metadata.json`. May leave header-only output on failure. |
| `-e`, `--essential-genes` | Required existing file | First line skipped; essential reference gene labels. |
| `-n`, `--non-essential-genes` | Required existing file | First line skipped; nonessential reference gene labels. |
| `-k`, `--use-column` | Exact header string; `BF` | Selects ranking score, e.g. `BoostedNormBF` for corrected network + guide normalization. Absent header fails. |
| `--dump-json` | Optional JSON output path | Writes raw PR rows and provenance. |
| `-h`, `--help` | Flag | Prints help. |

PR output columns are `Gene`, selected score name, `Recall`, `Precision`, `FDR`. Rows are sorted by descending selected score. Precision/FDR are reference-label estimates, not probabilities that each individual gene is essential.


To inspect the installed version:

```sh
./bin/rsbagel --version
```

Expected output: `rsbagel 0.1.0`. Success means the program exits normally after this line.

To list available commands:

```sh
./bin/rsbagel --help
```

Expected output: a usage line and the `fc`, `bf`, `pr`, `version`, and `help` commands. Success means the help text appears without an error.

To inspect fold-change options:

```sh
./bin/rsbagel fc --help
```

Expected output: usage for `rsbagel fc` and its options. Success means the help text includes `--read-count-file`, `--output-label`, and `--control-columns`.

To inspect scoring options:

```sh
./bin/rsbagel bf --help
```

Expected output: usage for `rsbagel bf` and its options. Success means the text includes the generator, worker, and network-profile choices.

To inspect ranking options:

```sh
./bin/rsbagel pr --help
```

Expected output: usage for `rsbagel pr` and its options. Success means the text includes `--use-column`.

Output folders must already exist. `fc -o` is a prefix to which suffixes are added. `bf -o` and `pr -o` are the exact main output file names. Existing output files can be replaced. Use a new folder for each analysis.

## Sampling profiles and experimental native mode

A random generator makes the sequence used to divide or resample genes. RSBagel offers two generators. “NumPy-compatible” describes sampling behavior; it does not require running Python.

| Choice | Status | Appropriate use |
|---|---|---|
| `--rng-compat numpy` | Production default within the documented qualified scope | Routine supported analyses; reproducing the pinned BAGEL2 sampling contract |
| `--rng-compat rust` | Experimental native ChaCha8; not production-qualified | Explicit exploratory or diagnostic work that is labelled experimental |

Use NumPy-compatible mode for supported compatibility work. The reference is BAGEL2 build115 at the pinned source revision. The default also applies to library integrations, described in the developer guide.

Native ChaCha8 is a valid deterministic implemented mode. It is not the default compatibility profile and is not an invalid or failed algorithm. It is available for explicit experiments. Its historical HAP1 comparison had Spearman correlation 0.999958347522601 and top-100 overlap 97/100. That overlap missed the required 98/100 threshold. High rank agreement does not make this mode production-qualified or biologically superior.

Optional experimental example:

```sh
./bin/rsbagel bf -i outputs/results.foldchange -o outputs/native.bf -e inputs/essentials.tsv -n inputs/nonessentials.tsv -c 1,2 --seed 42 --workers 1 --rng-compat rust
```

This explicitly selects native ChaCha8. Expected output begins `Using column:  R1, R2` and reports cross-validation iterations. Success creates `native.bf` and its metadata file. These files are experimental results. They are not expected to match the NumPy-compatible expected files.

The random-generator choice is separate from `--network-profile`. A network profile selects network-score behavior. It does not change the qualification of the random generator.

## Replay manifests

Replay is an advanced diagnostic input, not a routine requirement. The manifest
is JSON describing the exact gene order and train/test partitions. Use an
existing valid manifest for the exact input ordering and requested iteration
count. Do not treat replay as an interchangeable seed. Gene identity/order,
contiguous iteration IDs, complete disjoint membership, and iteration count
must agree. Bad JSON or inconsistent membership fails. The existing checked
fixture has 100 iterations; it requires a matching iteration request.

## Resampling modes and precedence

Choose one mode. Routine commands use default ten-fold cross-validation.
Bootstrap is selected with `--use-bootstrapping`; its default iteration count is
1000. Small-sample uses 100 reference-gene bootstrap fits. No-resampling trains
and scores all genes once and does not provide independent validation.

The parser permits combined flags. Small-sample takes precedence over bootstrap,
which takes precedence over cross-validation. No-resampling overrides partition
generation but retains the selected output schema. Avoid combined mode flags.

## Network options and risky legacy combinations

Networks are optional headerless tab-separated pairs of gene identifiers.
Selecting a profile alone does not enable a network. Gene-level network scores
require a network file; reagent-level output skips network correction.

**Advanced legacy behavior: not recommended for routine use.** Legacy mode keeps
three known interactions: no resampling uses a fresh network fold plan;
nondefault cross-validation counts interact with a fixed-ten-fold network plan;
and guide equalization mixes raw and normalized boosted observations. These can
produce nonfinite values, partial files, or failures. A zero exit status with
nonfinite scores is not a valid biological result. These defective combinations
have no scientific-validity claim.

Corrected mode reuses main partitions, records an individual missing network
contribution as zero only when the model is valid, and keeps `BoostedBF` and
`BoostedNormBF` separate. An invalid model still fails. Correction does not
establish biological superiority or A375 network validation.

The PR command defaults to `BF` even when BF metadata recommends another ranking
column. Select the actual desired score header explicitly. In corrected network
analysis with guide equalization, the recommended header is `BoostedNormBF`.

## Multi-target filtering and alignment files

Filtering requires `--filter-multi-target` and a matching `--align-info` file.
The alignment input has no header and exactly five tab-separated fields per row:
reagent ID, perfect-match loci, perfect-match genes, mismatch loci, mismatch
genes. Multiple entries in a field are comma-separated. IDs must be unique.
The default locus thresholds are 10 perfect and 10 mismatch loci. A missing
required alignment file or malformed alignment record fails. A supplied path
must exist even if filtering is off, because metadata records input identity.
Do not apply individual-guide mapping assumptions to component spacers by
splitting a measured array into artificial observations.

## Diagnostic outputs, raw JSON, and technical output fields

The routine example produces four scientific tables and three metadata files.

| File | Contents |
|---|---|
| `results.foldchange` | Reagent IDs, gene labels, and log2 fold changes; three decimal places |
| `results.normed_readcount` | Reagent IDs and normalized read-count columns; two decimal places; no gene-label column |
| `results.bf` | Gene scores, alphabetically ordered by gene label; three decimal places |
| `results.pr` | Genes ranked by the selected score, with reference-based precision and recall |
| `results.metadata.json` | Fold-change configuration, count-file identity, build and invocation provenance |
| `results.bf.metadata.json` | Scoring configuration, input identities, profile, ranking hint, and provenance |
| `results.pr.metadata.json` | Ranking score column, input identities, and provenance |

**JSON** is a structured text format used to store named settings and values. A metadata file is a record about how an analysis was made. Retain it even if you mainly inspect the tables.

Gene-level score fields depend on the requested options:

| Field | Meaning |
|---|---|
| `GENE` | Gene identifier |
| `BF` | Mean gene Bayes-factor score across scoring iterations; replicate equalization can scale it |
| `STD` | Standard deviation across bootstrap or small-sample gene-score iterations; not a biological-replicate error bar |
| `NumObs` | Number of scoring iterations contributing to the score; not the number of guides or sequenced cells |
| `NormBF` | Guide-count-equalized gene score when requested |
| `BoostedBF` | Score with network contribution when a network is used |
| `STD_BoostedBF` | Variation across boosted scoring iterations in bootstrap or small-sample output |
| `BoostedNormBF` | Separate network-boosted, guide-equalized score under the corrected profile |

Default cross-validation output has `GENE` and `BF`. Bootstrap and small-sample schemas add `STD` and `NumObs`. Optional network and equalization choices add their corresponding fields. Replicate equalization preserves legacy field-specific scaling; it does not multiply every reported field by the same factor.

The precision–recall table contains:

| Field | Meaning |
|---|---|
| `Gene` | Gene identifier |
| Selected score name | The score used for ranking; `BF` by default |
| `Recall` | Essential reference genes encountered so far divided by essential reference genes present in the input score table |
| `Precision` | Essential references encountered so far divided by essential plus nonessential references encountered so far |
| `FDR` | One minus the reported precision |

Unlabelled genes appear in the ranking but do not add to either reference counter. Before any reference is encountered, the implementation starts precision at 1. Scores tied after three-decimal serialization use pinned compatibility ordering, not an invented alphabetical tie rule.

Optional diagnostic outputs include `.traininfo` from test mode, `.netscore` from gene-level network test mode, and the user-named JSON file from `--dump-json`. Training files describe partitions. Network diagnostic rows contain gene, network score, and contribution. Raw JSON can contain `null` for nonfinite numeric values. A nonfinite value is not a valid zero score.

## Parser edge cases and numerical limits

Routine analyses should use complete nonnegative counts, a positive pseudocount,
positive column positions or exact names, and positive workers or `auto`.

Advanced compatibility behavior retained from validation:

- Reference files always lose their first line to the header skip. Headerless
  lists therefore silently lose the first gene.
- Count parsing fills supported missing numeric tokens with zero and missing
  gene tokens with `NO_GENE_NAME`. BF numeric parsing does not do that filling.
- Numeric-column indices change after controls are removed. Negative selectors
  have legacy stage-specific behavior; BF negative `-1` fails.
- Repeated value options keep the last value. `--compat` selects NumPy mode even
  beside an explicit native selection.
- Minimum reads and pseudocount use signed integers. Equalization accepts signed
  integers too; zero disables it, while negative values have no recommended
  scientific use. Fractional or noninteger strings are not accepted here.
- NumPy seeds range from 0 to 4294967295. Native seeds range from 0 to
  18446744073709551615. A NumPy seed above its range fails during computation.
- Worker values must be `auto` or positive integers. Explicit worker values have
  no automatic eight-worker cap.
- Zero CV folds or active bootstrap iterations fails. Excessive CV folds can
  fail for the input gene universe. Constant or insufficient training data can
  make numerical fitting fail. No universal minimum valid reference size is
  established.
- PR requires an existing score header and matching essential references. It can
  leave a header-only output on failure. Ties follow pinned ordering, not a
  stable alphabetical rule. Nonfinite scores are not zeros.

These are existing parser/domain contracts, not new qualification. Integer
boundaries are parser type limits; the earlier checks did not enumerate every
possible integer or malformed input. Exit 0 means command success, 1 means a
data/computation error, and 2 means a command syntax error. Existing output files
can be overwritten and failures can leave partial files. Keep a separate output
folder and never infer success merely from file existence.

## Performance, memory, and worker details

Start with the supplied small example. Use ordinary settings before adding diagnostic output, networks, or large bootstrap runs.

The default `--workers auto` uses the smaller of detected performance cores and eight on Apple Silicon. If topology is unavailable, it uses the smaller of available parallelism and eight, with a minimum of one. A positive worker count overrides this. A worker is one unit of parallel work. `--workers 1` is useful for diagnosis and is used in the teaching example.

More workers are not a guarantee of faster completion for small inputs. Bootstrap iteration count, reagent count, selected sample count, and network processing affect time and memory. `--dump-json` and test mode retain additional internal information. Use them only when needed. Close other memory-heavy applications before a large job and retain enough disk space for outputs.

Recorded performance measurements apply to their stated workloads and hardware. They do not establish universal speedups. Memory use is also workload-dependent. Real-data cases showed reductions in the documented comparisons, while the large B4 synthetic workload used about 1.56–1.59 times the Python peak memory. Do not assume RSBagel always uses less memory. This manual does not prescribe a universal RAM minimum or promise performance on unqualified machines.

## Advanced troubleshooting

| Symptom | Likely cause | Action |
|---|---|---|
| `No such file or directory` after `cd` | Package is not at the expected Desktop location | Move the complete folder to Desktop in Finder; keep the exact `RSBagel-example` name |
| Program or script not found | Wrong working directory or incomplete folder | Repeat the folder-selection step; confirm `bin/rsbagel` and `run-example.sh` are present |
| `Permission denied` | Executable permissions were lost during copying | Ask the package provider to restore executable permissions or provide the intact package |
| macOS blocks execution | Local security policy or a quarantined, unsigned file | Ask your administrator or package provider to verify this local binary; do not disable system security globally |
| Wrong architecture or `Bad CPU type` | Wrong executable for this Mac | Use the Apple Silicon executable on the documented platform |
| `Invalid input controls` or `Invalid columns` | Incorrect selector | Count numeric columns from 1; recheck column names after `fc` removes controls |
| Could not convert text to float | Text in a numeric cell or the wrong input format | Re-export tab-separated text and inspect that row |
| Unequal field counts | Commas, embedded tabs, or missing cells | Check delimiter and column count in every row |
| Missing control genes or numerical model failure | Gene names do not match, too little control data, or degenerate distributions | Check lists, headers, labels, and control coverage; do not interpret failed fits |
| Multi-target filtering requests alignment information | Filter enabled without its required alignment file | Supply a valid matching file or omit the filter |
| PR says the column name is absent | Selected score is not in that BF table | Inspect the header and use an actual score field |
| PR reports division by zero | No essential reference genes match the score table | Fix gene labels and reference list; a partial output header is not success |
| Example expected-file comparison fails | Changed inputs, program, settings, or incomplete run | Restore the supplied package and repeat the example; retain the failed files for diagnosis |
| Different results after rerunning | Omitted seed, changed generator, files, or settings | Restore the full recorded inputs and configuration, including the generator |
| `nan`, `inf`, or empty numeric cells | Nonfinite model result or known legacy combination | Stop interpretation and diagnose the model and profile; never substitute zero silently |

Check the error message as well as the output folder. `bf` or `pr` can leave partial diagnostic files after a failure. Treat a nonzero exit as a failed run, even when some files exist.

## Detailed reproducibility contract

A **seed** is the integer that initializes random sampling. Give an explicit seed when you want to reproduce a run. The CLI generates a time-derived seed when `--seed` is omitted. The Rust core's default seed is 42; do not confuse that library default with the CLI's behavior.

Equal numerical seeds do **not** produce equal random sequences under different generators. Seed 42 with NumPy-compatible sampling and seed 42 with ChaCha8 represent different sampling sequences.

Keep these items together:

- The executable version and identity.
- All input files, including reference lists.
- The complete command and selected columns.
- The seed and generator.
- The resampling settings.
- The network profile and any optional graph or alignment file.
- The metadata sidecars.

Scientific reductions use a fixed order. Worker-count reproducibility was tested within the documented platform and cases. This is not a guarantee of universal internal bitwise identity on every platform, compiler, input, or operation. The demonstrated A375 compatibility concerns serialized scientific outputs and ordinal rankings on those datasets.

The teaching package compares scientific text files byte for byte. Metadata files contain paths, command arguments, configuration, and build provenance. Their bytes can change when the folder location or invocation changes. Do not require metadata bytes to remain identical across different folders.

## Version, qualification limits, citation, and notices

This manual covers RSBagel 0.1.0 and the current inspected source. The `version` subcommand also prints the BAGEL2 compatibility target, Version 2.0 / Build 115. That target is not the RSBagel software version.

Compatibility evidence is bounded by the documented synthetic, HAP1, and A375 cases and their settings. Exact A375 scientific text and rank agreement does not mean every internal floating-point operation is universally bitwise identical. Native ChaCha8 remains experimental and unqualified. No new biological discovery is claimed.

The legacy network profile intentionally retains three known interactions: no resampling with a fresh network fold plan, nondefault cross-validation counts with the fixed-ten-fold network plan, and mixed raw/normalized boosted observations under guide equalization. These defective combinations have no scientific-validity claim. The corrected network profile reuses main partitions, records missing per-gene evidence as a zero contribution only when the model is valid, and separates boosted raw and normalized scores. Invalid network models still fail. Correction does not establish biological superiority or A375 network validation.

For the underlying method, the accepted bibliography lists Hart and Moffat (2016), BAGEL, DOI 10.1186/s12859-016-1015-8, and Kim and Hart (2021), BAGEL2, DOI 10.1186/s13073-020-00809-3. The accepted RSBagel manuscript title is “RSBagel Preserves the BAGEL2 Scientific Contract while Accelerating Pooled CRISPR Fitness-Screen Analysis.” Treat it as a local manuscript, not an already published article. Record the RSBagel version, source identity, and selected profile with your methods. A public RSBagel release DOI and finalized publication citation are not established by this local package. Do not invent one. Acknowledge the underlying BAGEL2 method using the accepted manuscript's BAGEL2 reference.

RSBagel uses the MIT license. BAGEL2-derived behavior carries the Hart Lab MIT notice. Pinned NumPy indirect-sort semantics carry a BSD notice. The package includes the corresponding license texts. Keep these notices with redistributed copies.

## Reagent-output schema details

In cross-validation reagent output, columns start with `RNA`, `GENE`, selected
sample names, and `BF`. `RNA` is the reagent identifier, not a gene name.
Bootstrap and small-sample output also contain sample `_STD` columns and
`NumObs`. Keep array rows whole in this output too. This table is for inspecting
reagents and is not a gene-level PR input.
