# RSBagel Quick Start

Use an Apple Silicon Mac and an RSBagel v1.0.0 source checkout. Open Terminal and run each command below from the repository root.

Use the default build-output location:

```sh
unset CARGO_TARGET_DIR
```

Build the program:

```sh
cargo build --release --locked --target aarch64-apple-darwin
```

Create the example’s `bin` folder:

```sh
mkdir -p examples/beginner/bin
```

Copy the built executable into the example:

```sh
cp target/aarch64-apple-darwin/release/rsbagel examples/beginner/bin/rsbagel
```

## Copy the folder

In Finder, copy `examples/beginner` to your `Desktop`. Name the copied folder exactly `RSBagel-example`. Keep Terminal open for the next steps.

## Run the example

1. Paste this command into Terminal and press Return.

```sh
cd "$HOME/Desktop/RSBagel-example"
```

This selects the example folder. `$HOME` means your home folder; keep it as
written. Expected output: none. Success means the prompt returns without an
error.

2. Paste this command and press Return.

```sh
./run-example.sh
```

This runs the example and checks its results. It prints progress messages.
Success means the last two lines are:

```text
SUCCESS: all four scientific result files match expected outputs.
Results are in outputs/results.pr and outputs/results.bf.
```

## Find the results

1. Open `RSBagel-example` on your `Desktop` in Finder.
2. Open its `outputs` folder.
3. Open `results.pr` in a text editor or import it into a spreadsheet using Tab as the separator.
4. Retain all three files whose names end in `.metadata.json`.

`results.pr` is the ranked result table. `results.bf` contains gene scores.
For explanations and your own data, use the [User Manual](RSBagel-User-Manual.md).

## Fix common first-run errors

| Message or problem | Fix |
|---|---|
| `no such file or directory` | Put the complete folder on `Desktop` and check its exact name. Repeat the folder-selection command. |
| `permission denied` | Rebuild and copy the executable as shown above; check that the script and executable have execute permission. |
| macOS blocks the program | Follow your institution's approved process for local software. |
| `Bad CPU type in executable` | Use an Apple Silicon Mac. |
| `outputs already exists` | Move the old `outputs` folder elsewhere in Finder before repeating the run. |
| Missing input or no `SUCCESS` message | Keep the error text. Restore the intact package and retry. Do not use incomplete results. |
