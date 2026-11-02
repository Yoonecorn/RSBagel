# RSBagel User Manual

RSBagel 1.0.0. Start with the [Quick Start](RSBagel-Quick-Start.md) to run the supplied example.

Most users need only three steps:

1. Prepare the count table and two control-gene lists.
2. Run `fc`, `bf`, and `pr`.
3. Review `results.pr` and retain all metadata files.

Go directly to [Use your own data](#use-your-own-data) once the example succeeds.

## Overview

RSBagel analyzes pooled CRISPR fitness screens. It compares changes in reagent
counts with essential and nonessential control genes. It reports gene scores
and a ranked result table. Higher scores support essential-like depletion in
that screen. They do not establish that a gene is essential in every setting.

A command-line interface (CLI) lets you run the program by typing commands in
Terminal. The executable, sometimes called a binary, is the program file
`bin/rsbagel`. A working directory is the folder Terminal uses for file names.
Here it is `RSBagel-example` on your `Desktop`.

## Files you need

- [ ] Count table
- [ ] Essential-gene list
- [ ] Nonessential-gene list
- [ ] RSBagel executable
- [ ] New output folder

Network, alignment, replay, and JSON files are not required as inputs for a
routine analysis. RSBagel writes metadata files automatically; retain them.
Metadata means the recorded settings and input identities for a run.

## Installation

This package runs on an Apple Silicon Mac. First build the executable and prepare the example folder using the [Quick Start](RSBagel-Quick-Start.md). Python and NumPy are not needed to run the built executable.

1. Copy the complete supplied `RSBagel-example` folder to your `Desktop` in Finder.
2. Keep its folder name unchanged.
3. Press Command–Space.
4. Type `Terminal`.
5. Press Return.

Terminal is the app where you type commands. Select the example folder:

```sh
cd "$HOME/Desktop/RSBagel-example"
```

`$HOME` means your home folder; leave it as written. It is an environment variable,
a named setting available to commands. This command selects the working
directory. Expected output: none. Success means the prompt returns without an
error. The explicit program path `./bin/rsbagel` needs no PATH setup. PATH is the
list of folders searched for a program name. `./` means the current folder.

Source qualification is limited to the evaluated Apple Silicon macOS configuration. See the [manuscript](RSBagel-Manuscript.pdf) for the validation environment and support boundaries.

## First analysis

### Run the supplied example

```sh
./run-example.sh
```

This runs all three analysis stages and compares the scientific tables with the
supplied expected files. Progress messages are normal. Success ends with:

```text
SUCCESS: all four scientific result files match expected outputs.
Results are in outputs/results.pr and outputs/results.bf.
```

Open `outputs/results.pr` in a text editor or import it as a tab-separated table
in a spreadsheet. Retain all three `.metadata.json` files in `outputs`.
The example contains artificial data: 120 genes, three separately measured
reagents per gene, a starting sample, and two endpoint samples.

### Run the same stages separately (optional)

These commands use the example files. Run them after the example has created
`outputs`. They replace its existing result files. For your data, use the new
output folder described in the next workflow.

1. Calculate fold changes.

```sh
./bin/rsbagel fc -i inputs/counts.tsv -o outputs/results -c 1
```

This selects the starting count column `T0`, the first numeric column.
Expected output:

```text
Controls: T0
```

Success creates `results.foldchange`,
`results.normed_readcount`, and `results.metadata.json` in `outputs`.

2. Calculate gene scores.

```sh
./bin/rsbagel bf -i outputs/results.foldchange -o outputs/results.bf -e inputs/essentials.tsv -n inputs/nonessentials.tsv -c 1,2 --seed 42 --workers 1
```

This scores both endpoint columns with the two control lists. A seed starts
random sampling; seed 42 makes repeat sampling possible. A worker is one unit
of parallel calculation; this command uses one. The beginner workflow uses one
worker to provide a simple and reproducible example. Expected output starts with:

```text
Using column:  R1, R2
```

It includes ten scoring iterations. Success creates
`results.bf` and `results.bf.metadata.json`.

3. Create the ranked table.

```sh
./bin/rsbagel pr -i outputs/results.bf -o outputs/results.pr -e inputs/essentials.tsv -n inputs/nonessentials.tsv
```

This ranks genes by their Bayes-factor score (BF). Expected Terminal output:
none. Success creates `results.pr` and `results.pr.metadata.json`.

### Use your own data

This copy-and-paste workflow assumes **one starting-sample column and two
endpoint columns**, like the example. For another number of endpoints, replace
`1,2` with the positive positions of the endpoint columns you want to score in
the fold-change file. The starting column has been removed from that file.

The capitalized names below are placeholders: replace them with your own file,
column, and folder names before pasting. Keep the quotes. The definitions follow
each command. Run one command at a time and stop if it reports an error.

1. Copy your count table and both control-gene lists into the `RSBagel-example` folder in Finder.
2. Create a new, empty output folder there in Finder.
3. Replace the placeholders in the fold-change command.
4. Run the fold-change command.

```sh
./bin/rsbagel fc -i "MY_COUNTS.tsv" -o "MY_ANALYSIS/results" -c "MY_STARTING_SAMPLE"
```

- `MY_COUNTS.tsv`: your tab-separated count-table file name in the working folder.
- `MY_ANALYSIS`: the new empty output folder you just created; keep `/results` as written.
- `MY_STARTING_SAMPLE`: the exact header of the starting-sample count column, such as `T0`.

This calculates fold changes with the default positive pseudocount of 5.
Expected output starts with this prefix, followed by your starting-column name:

```text
Controls:
```

Success creates `results.foldchange`, `results.normed_readcount`, and
`results.metadata.json` in your chosen folder.

5. Replace the placeholders in the scoring command.
6. Run the scoring command.

```sh
./bin/rsbagel bf -i "MY_ANALYSIS/results.foldchange" -o "MY_ANALYSIS/results.bf" -e "MY_ESSENTIALS.tsv" -n "MY_NONESSENTIALS.tsv" -c 1,2 --seed 42 --workers 1
```

- `MY_ANALYSIS`: the same output folder used in the previous command.
- `MY_ESSENTIALS.tsv`: your essential-gene list file name in the working folder.
- `MY_NONESSENTIALS.tsv`: your nonessential-gene list file name in the working folder.

This uses the first two numeric fold-change columns, seed 42, one worker,
the default NumPy-compatible profile, and only the default ten-fold
cross-validation resampling mode. Cross-validation divides genes into groups
so each group can be scored using a model fitted to other groups.
Expected output starts with this prefix and your two endpoint names:

```text
Using column:
```

Ten scoring iterations follow. Success creates `results.bf` and
`results.bf.metadata.json`. Do not add other resampling flags to this routine workflow.

7. Replace the placeholders in the ranking command.
8. Run the ranking command.

```sh
./bin/rsbagel pr -i "MY_ANALYSIS/results.bf" -o "MY_ANALYSIS/results.pr" -e "MY_ESSENTIALS.tsv" -n "MY_NONESSENTIALS.tsv"
```

- `MY_ANALYSIS`: the same output folder used in both earlier commands.
- `MY_ESSENTIALS.tsv`: the same essential-gene list used for scoring.
- `MY_NONESSENTIALS.tsv`: the same nonessential-gene list used for scoring.

This ranks gene scores and reports precision and recall. Expected Terminal
output: none. Success creates `results.pr` and `results.pr.metadata.json`.
Open your chosen output folder to find the results. Keep every metadata file.

## Input preparation

Use plain text with tabs between columns. The extension `.tsv` means tab-separated values. Do not provide an Excel workbook directly.

A count table begins with a header row:

```tsv
REAGENT_ID	GENE	T0	R1	R2
geneA_sg1	geneA	1000	100	120
geneA_sg2	geneA	900	130	110
geneB_sg1	geneB	1100	1000	1030
```

This table is a format illustration, not a complete training dataset.

1. Give every observed reagent a unique identifier.
2. Put the target gene in the second column.
3. Put the starting-sample counts in a numeric column.
4. Put each later sample in its own numeric column.
5. Use distinct sample names.
6. Export as tab-separated plain text.
7. Check that every row has the same number of fields as the header.

Use complete, nonnegative integer read counts. Supply actual measured counts. Do not replace missing samples with fabricated values. Keep identifiers and gene labels consistent across files. Extra spaces in labels can prevent matching.

`fc -c` selects starting-sample columns. Numeric selection counts only the numeric columns: `T0` is `1`, `R1` is `2`, and `R2` is `3` in the illustrated count table. The gene column is not a count column. A comma-separated list selects more than one control column; those control counts are combined. Exact column names can also be used. Avoid negative selectors and repeated selectors in routine work.

The output fold-change table contains `REAGENT_ID`, `GENE`, and the remaining numeric sample columns. Control columns are removed. Consequently, `R1` becomes column `1` and `R2` becomes column `2` for `bf -c` in the example. Reusing the original count-table indices at this stage is a common error.

### Technical note: how RSBagel calculates fold change

Normalization first sums selected control counts, then adds the pseudocount once to that combined control and once to each remaining sample. It rescales each resulting column to 10 million counts. The default pseudocount is 5. Fold change is log2 of normalized sample abundance divided by normalized control abundance. The minimum-read filter compares combined control abundance **after adding the pseudocount** with the threshold. The default threshold is zero. Preserve these settings when reproducing a previous run.

## Control-gene preparation

Two different kinds of controls are needed. The starting-sample columns define the abundance baseline. The essential and nonessential gene lists train the scoring model. These are not interchangeable.

Each gene list is a plain-text file with a header on the first line and one gene per later line. The parser skips the first line. A missing header silently loses the first gene.

```text
GENE
geneA
geneC
geneD
```

1. Choose reference lists appropriate for the organism and screen design.
2. Match the gene identifiers exactly to the count table.
3. Include a header in each reference file.
4. Keep the essential and nonessential lists separate.
5. Check that enough genes from both lists occur in the data.
6. Retain the exact lists with the analysis.

The example uses artificial control labels. Do not use them for human screen analysis. A non-targeting reagent label is not automatically a nonessential-gene reference. Small or constant control distributions can make model fitting fail. No universal minimum reference-list size is established here; the teaching package is known to work with its supplied inputs.

The references also determine the precision–recall estimates. Using them to evaluate a ranking does not provide an independent biological validation set.

## Routine supported profile

Use the default NumPy-compatible profile unless an advanced analysis specifically requires the experimental native profile.

The default needs no NumPy installation.

Native ChaCha8 is an implemented deterministic mode. It is not the default compatibility profile and remains experimental.

For executable options, use `./bin/rsbagel bf --help`. Equal numerical seeds do not make the native and NumPy-compatible random streams identical.

Use the positive pseudocount of 5, seed 42, default ten-fold cross-validation,
and a positive worker count or `auto` for routine work. Use positive column
positions or exact names. The commands above use one worker. Do not add network,
alignment, replay, or diagnostic options to complete this routine workflow.

## Reproducibility

Keep the executable, input files, complete command, seed, and metadata together.
Use an explicit seed such as 42. Without one, the CLI chooses a time-derived seed.
Equal numerical seeds do not give equal random sequences with different generators.
Keep the same generator, selected columns, and resampling settings when repeating
an analysis.

The example compares scientific result files byte for byte. Metadata includes
paths and command arguments, so those files may differ after moving the folder.
There is no universal promise of internal bitwise identity across platforms or
compilers. Bitwise identity means identical stored numerical bits.

## Output files

The routine workflow creates four scientific tables and three metadata files.
All are inside `outputs` for the example, or your chosen output folder for your data.

| File | What it contains |
|---|---|
| `results.foldchange` | Each reagent's log2 abundance changes and gene label |
| `results.normed_readcount` | Normalized counts; reagent IDs but no gene-label column |
| `results.bf` | Gene scores in alphabetical gene order |
| `results.pr` | Genes ranked from larger to smaller BF, with recall, precision, and FDR |
| `results.metadata.json` | Starting-sample settings and count-file identity |
| `results.bf.metadata.json` | Scoring settings and input identities |
| `results.pr.metadata.json` | Ranking settings and input identities |

JSON is a structured text format for named values. You do not need to prepare
JSON inputs. Retain the automatically written metadata files with the results.
For optional command-line diagnostics, consult the relevant command’s `--help` output.

## Result interpretation

A negative log2 fold change means relative depletion. A value near −1 means about half the normalized abundance of the starting sample. Positive values mean relative enrichment. Library normalization means these are relative abundance comparisons, not direct cell counts.

A larger BF supports the essential-like model over the nonessential-like model. It is not a probability, a p-value, or a fold change. A small or negative BF does not prove that a gene has no biological role. Read scores with screen quality, reference behavior, guide agreement, and experimental design in mind.

The PR table is sorted by decreasing score. Choose a score threshold only after considering reference performance and the purpose of the experiment. Its FDR is reference-based. It is not a universal guarantee of the proportion of false discoveries among all unlabelled genes.

Do not compare raw BF values across screens with different reagent numbers or replicate selections as though their scales were automatically identical. Optional equalization changes scale and interpretation. Do not claim that one guide system and another are biologically equivalent because their software outputs agree.

### Worked example: the first ranked row

The unchanged validated example begins with this row:

| Field | Value | Meaning |
|---|---:|---|
| Gene | G010 | An artificial gene in the simulated essential-control list |
| BF | 878.089 | Strong support for an essential-like depletion pattern relative to the nonessential-like model; not a probability or p-value |
| Recall | 0.025 | One of the 40 measured essential controls has been reached: 1/40 = 2.5% |
| Precision | 1.000 | All reference controls encountered so far are essential: one essential and no nonessential controls |
| FDR | 0.000 | One minus precision at this point; no nonessential reference has yet been encountered |

This row supports a high essential-like score and a favorable first-ranked
reference-control result in this synthetic example. It does not prove a
biological effect. It does not mean an unlabelled gene with the same score would
have zero probability of being a false discovery.

## Individual-sgRNA workflow

A Brunello-style individual sgRNA is a separate observed reagent. Several sgRNAs targeting the same gene occupy separate rows with different reagent IDs and the same gene label. Keep their measured counts separate. The routine workflow combines their evidence at gene level.

Optional reagent-level output helps inspect agreement between observed sgRNAs:

```sh
./bin/rsbagel bf -i outputs/results.foldchange -o outputs/guides.bf -e inputs/essentials.tsv -n inputs/nonessentials.tsv -c 1,2 --seed 42 --workers 1 -r
```

This scores the same data and writes reagent-level rows. Expected output begins:

```text
Using column:  R1, R2
```

The iteration table follows. Success creates `guides.bf` and `guides.bf.metadata.json`. The table begins `RNA`, `GENE`, the selected sample names, and `BF` for this cross-validation example. `RNA` is the reagent identifier, not a gene name. Additional output fields are explained in the advanced reference.

Keep this table for guide inspection. Do not pass it directly to the gene-level PR workflow: the first column contains reagent IDs, and repeated gene labels do not turn reagent scores into a gene score table. Run gene-level `bf` separately to obtain the gene-level ranking.

## Four-guide single-target array workflow

An Inzolia-style four-guide single-target array is **one observation**. A count measures the whole observed array, not each component guide independently.

1. Keep one row for each measured array.
2. Use the array identifier as `REAGENT_ID`.
3. Use the single admitted target gene as `GENE`.
4. Preserve the measured array counts unchanged.
5. Use the ordinary gene-level workflow with those rows.

Never divide one array into four artificial guide observations. Do not copy the same count into four rows. Do not divide the count by four. Both would invent observations that the experiment did not measure.

If several separately measured single-target arrays target one gene, retain each observed array as its own row. An array that targets several different genes does not meet this single-target workflow. Resolve its interpretation before preparing a gene-level input. The included teaching input uses individual artificial reagents; it is not a second validated array dataset.

## Routine troubleshooting

| Problem | Action |
|---|---|
| The folder-selection command fails | Check that the complete `RSBagel-example` folder is on `Desktop`. |
| Program not found or permission denied | Check the working folder and ask the provider for an intact executable copy. |
| macOS blocks the program or reports a wrong CPU type | Use an Apple Silicon Mac and follow your institution's local-software process. |
| The example says `outputs already exists` | Move the previous `outputs` folder elsewhere in Finder before repeating it. |
| A real-data output file already exists | Use a new empty output folder to avoid replacing earlier results. |
| An input file cannot be found | Check its file name and location; replace every placeholder before running. |
| Invalid controls or columns | Match the starting header exactly; count endpoint columns from 1 after the starting column is removed. |
| A number cannot be read or row widths differ | Export plain tab-separated text and check the affected row. |
| Training fails or PR reports division by zero | Check headers, matching gene labels, and usable reference genes from both lists. |
| No success message, incomplete files, or nonfinite numbers | Stop interpretation and keep the error text. Restore the intact example or ask for help checking the data. |
| Results change after a repeat | Check the exact files, generator, seed, columns, and settings. |

A file can remain after an error. Stop if a command reports an error; file
existence alone is not success. For specialist errors, retain the exact command and error text with the run metadata.

## Limitations, citation, and license

This manual covers RSBagel 1.0.0.
Compatibility evidence applies to the documented inputs and settings. It does
not prove biological truth or universal internal numerical identity. Native
ChaCha8 retains its experimental qualification boundary. Performance and memory
use depend on the workload; no universal memory reduction or RAM minimum is
promised. Additional technical details are in the
[Developer Integration Guide](../developer/RSBagel-Developer-Integration-Guide.md).

For the underlying method, the accepted bibliography lists Hart and Moffat (2016), BAGEL, DOI 10.1186/s12859-016-1015-8, and Kim and Hart (2021), BAGEL2, DOI 10.1186/s13073-020-00809-3. The accepted RSBagel manuscript title is “RSBagel Preserves the BAGEL2 Scientific Contract while Accelerating Pooled CRISPR Fitness-Screen Analysis.” Treat it as a local manuscript, not an already published article. Record the RSBagel version, source identity, and selected profile with your methods. A public RSBagel release DOI and finalized publication citation are not established by this local package. Do not invent one. Acknowledge the underlying BAGEL2 method using the accepted manuscript's BAGEL2 reference.

RSBagel uses the MIT license. BAGEL2-derived behavior carries the Hart Lab MIT notice. Pinned NumPy indirect-sort semantics carry a BSD notice. The package includes the corresponding license texts. Keep these notices with redistributed copies.

Wet-lab users can skip the [Developer Integration Guide](../developer/RSBagel-Developer-Integration-Guide.md).
