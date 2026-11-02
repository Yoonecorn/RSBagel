use clap::{Args, Parser, Subcommand};
use rsbagel_core::{self as core, input, output, *};
use std::{
    fs::File,
    io::BufWriter,
    path::{Path, PathBuf},
};

#[derive(Parser)]
#[command(
    name = "rsbagel",
    version,
    about = "BAGEL2 f64 Rust compatibility engine"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}
#[derive(Subcommand)]
enum Command {
    Fc(FcArgs),
    Bf(Box<BfArgs>),
    Pr(PrArgs),
    Version,
}
#[derive(Args)]
#[command(args_override_self = true)]
struct FcArgs {
    #[arg(short = 'i', long = "read-count-file")]
    input: PathBuf,
    #[arg(short = 'o', long = "output-label")]
    output: PathBuf,
    #[arg(short = 'c', long = "control-columns", allow_hyphen_values = true)]
    controls: String,
    #[arg(short = 'm', long = "min-reads", default_value_t = 0)]
    min_reads: i64,
    #[arg(long = "pseudo-count", default_value_t = 5, allow_hyphen_values = true)]
    pseudo_count: i64,
    #[arg(long)]
    dump_json: Option<PathBuf>,
}
#[derive(Args)]
#[command(args_override_self = true)]
struct BfArgs {
    #[arg(short = 'i', long = "fold-change")]
    input: PathBuf,
    #[arg(short = 'o', long = "output-file")]
    output: PathBuf,
    #[arg(short = 'e', long = "essential-genes")]
    essential: PathBuf,
    #[arg(short = 'n', long = "non-essential-genes")]
    nonessential: PathBuf,
    #[arg(short = 'c', long = "columns-to-test", allow_hyphen_values = true)]
    columns: String,
    #[arg(short = 'w', long = "network-file")]
    network: Option<PathBuf>,
    #[arg(short = 'm', long = "filter-multi-target")]
    multi: bool,
    #[arg(long = "loci-without-mismatch", default_value_t = 10)]
    max_perfect: usize,
    #[arg(long = "loci-with-mismatch", default_value_t = 10)]
    max_mismatch: usize,
    #[arg(long = "align-info")]
    alignment: Option<PathBuf>,
    #[arg(short = 'b', long = "use-bootstrapping")]
    bootstrap: bool,
    #[arg(long = "no-resampling")]
    no_resampling: bool,
    #[arg(long = "use-small-sample")]
    small: bool,
    #[arg(short = 'N', long = "no-of-cross-validations", default_value_t = 10)]
    folds: usize,
    #[arg(long = "bootstrap-iterations", default_value_t = 1000)]
    iterations: usize,
    #[arg(short = 'r', long = "sgrna-bayes-factors")]
    guide: bool,
    #[arg(short = 'f', long = "equalise-sgrna-no", allow_hyphen_values = true)]
    equalise_guides: Option<i64>,
    #[arg(short = 'p', long = "equalise-rep-no", allow_hyphen_values = true)]
    equalise_reps: Option<i64>,
    #[arg(short = 's', long)]
    seed: Option<u64>,
    #[arg(short = 't', long = "run-test-mode")]
    test: bool,
    /// Sampling RNG: numpy is the production default; rust is experimental native ChaCha8,
    /// not production-qualified and excluded from compatibility speedup claims.
    #[arg(long="rng-compat",default_value="numpy",value_parser=["rust","numpy"])]
    rng: String,
    #[arg(long="network-profile",default_value="legacy",value_parser=["legacy","corrected"])]
    profile: String,
    #[arg(long = "replay-manifest")]
    replay: Option<PathBuf>,
    /// auto: performance cores (fallback: available parallelism), capped at 8.
    #[arg(long = "workers", alias = "threads", default_value = "auto", value_parser = core::workers::parse_workers)]
    threads: usize,
    #[arg(short = 'v', long)]
    verbose: bool,
    /// Explicit alias for --rng-compat numpy.
    #[arg(long)]
    compat: bool,
    #[arg(long)]
    dump_json: Option<PathBuf>,
}
#[derive(Args)]
#[command(args_override_self = true)]
struct PrArgs {
    #[arg(short = 'i', long = "bayes-factors")]
    input: PathBuf,
    #[arg(short = 'o', long = "output-file")]
    output: PathBuf,
    #[arg(short = 'e', long = "essential-genes")]
    essential: PathBuf,
    #[arg(short = 'n', long = "non-essential-genes")]
    nonessential: PathBuf,
    #[arg(short = 'k', long = "use-column")]
    column: Option<String>,
    #[arg(long)]
    dump_json: Option<PathBuf>,
}
fn writer(path: &Path) -> Result<BufWriter<File>> {
    Ok(BufWriter::new(File::create(path)?))
}
fn suffix(path: &Path, s: &str) -> PathBuf {
    let mut name = path.as_os_str().to_os_string();
    name.push(s);
    name.into()
}
fn dump(path: &Option<PathBuf>, value: &serde_json::Value) -> Result<()> {
    if let Some(path) = path {
        serde_json::to_writer(writer(path)?, value)?;
    }
    Ok(())
}
fn identity(path: &Path) -> Result<serde_json::Value> {
    use sha2::{Digest, Sha256};
    let bytes = std::fs::read(path)?;
    Ok(
        serde_json::json!({"path":path,"bytes":bytes.len(),"sha256":format!("{:x}",Sha256::digest(&bytes))}),
    )
}
fn provenance() -> serde_json::Value {
    serde_json::json!({"engine":"rsbagel","source_commit":core::SOURCE_COMMIT,"source_dirty":core::SOURCE_DIRTY,"version":env!("CARGO_PKG_VERSION"),"upstream_commit":"53388adbb4fb0931e5c9dda135502be19e4555f0","precision":"f64","reduction":"ordered iterations","argv":std::env::args().collect::<Vec<_>>()})
}
fn run(command: Command) -> Result<()> {
    match command {
        Command::Version => {
            println!(
                "BAGEL2 compatibility\nVersion: 2.0\nBuild: 115\nRSBagel {}",
                env!("CARGO_PKG_VERSION")
            );
        }
        Command::Fc(args) => {
            let model = input::read_matrix(File::open(&args.input)?, true)?;
            let (controls, labels) = input::select_columns(&model, &args.controls, false)?;
            println!("Controls: {}", labels.join(", "));
            let cfg = FcConfig {
                controls,
                control_labels: labels,
                min_reads: args.min_reads as f64,
                pseudocount: args.pseudo_count as f64,
            };
            let result = run_fc(&model, &cfg)?;
            output::write_matrix(
                writer(&suffix(&args.output, ".foldchange"))?,
                &result.fold_change,
                false,
            )?;
            output::write_matrix(
                writer(&suffix(&args.output, ".normed_readcount"))?,
                &result.normalized,
                true,
            )?;
            dump(
                &args.dump_json,
                &serde_json::json!({"config":cfg,"input":model,"result":result,"inputs":{"matrix":identity(&args.input)?},"provenance":provenance()}),
            )?;
            serde_json::to_writer_pretty(
                writer(&suffix(&args.output, ".metadata.json"))?,
                &serde_json::json!({"workflow":"fc","config":cfg,"inputs":{"matrix":identity(&args.input)?},"provenance":provenance()}),
            )?;
        }
        Command::Bf(args) => {
            if args.verbose {
                eprintln!("Resolved workers: {}", args.threads);
            }
            let model = input::read_matrix(File::open(&args.input)?, false)?;
            let ess = input::read_references(File::open(&args.essential)?)?;
            let non = input::read_references(File::open(&args.nonessential)?)?;
            let (columns, labels) = input::select_columns(&model, &args.columns, true)?;
            if args.multi && args.alignment.is_none() {
                return Err(Error::Invalid("Error! Multi-target-filtering selected and not align-info provided.\nPlease indicate align-info file.".into()));
            }
            let cfg = BfConfig {
                columns,
                column_labels: labels,
                resampling: if args.small {
                    Resampling::SmallSample
                } else if args.bootstrap {
                    Resampling::Bootstrap
                } else {
                    Resampling::CrossValidation
                },
                no_resampling: args.no_resampling,
                folds: args.folds,
                bootstrap_iterations: args.iterations,
                seed: args.seed.unwrap_or_else(|| {
                    ((std::time::SystemTime::now()
                        .duration_since(std::time::UNIX_EPOCH)
                        .unwrap_or_default()
                        .as_secs_f64()
                        * 100000.0)
                        % 100000.0) as u64
                }),
                rng_mode: if args.compat || args.rng == "numpy" {
                    RngMode::Numpy
                } else {
                    RngMode::Rust
                },
                threads: args.threads,
                replay: args
                    .replay
                    .as_ref()
                    .map(|p| -> Result<ReplayManifest> {
                        Ok(serde_json::from_reader(File::open(p)?)?)
                    })
                    .transpose()?,
                guide_level: args.guide,
                equalise_guides: args.equalise_guides.map(|v| v as f64),
                equalise_reps: args.equalise_reps.map(|v| v as f64),
                alignment: if args.multi {
                    Some(input::read_alignment(File::open(
                        args.alignment.as_ref().unwrap(),
                    )?)?)
                } else {
                    None
                },
                max_perfect_loci: args.max_perfect,
                max_mismatch_loci: args.max_mismatch,
                network: args
                    .network
                    .as_ref()
                    .map(|p| -> Result<_> { input::read_network(File::open(p)?) })
                    .transpose()?,
                network_profile: if args.profile == "corrected" {
                    NetworkProfile::Corrected
                } else {
                    NetworkProfile::Legacy
                },
                dump_stages: args.dump_json.is_some() || args.test,
            };
            if let Some(alignment) = &cfg.alignment {
                let discarded = alignment
                    .iter()
                    .filter(|a| {
                        a.perfect_loci > cfg.max_perfect_loci
                            || a.mismatch_loci > cfg.max_mismatch_loci
                    })
                    .count();
                println!("Total {discarded} multi-targeting gRNAs are discarded");
            }
            println!("Using column:  {}", cfg.column_labels.join(", "));
            let mut result = core::bf::run_bf_base(&model, &ess, &non, &cfg, None)?;
            if args.test {
                output::write_traininfo(writer(&suffix(&args.output, ".traininfo"))?, &result)?;
            }
            core::multitarget::apply_multitarget(&mut result, &cfg)?;
            let network_status = core::network::apply_network(&mut result, &ess, &non, &cfg);
            // Emit observed main-stage diagnostics before network-stage text; network edge count is
            // loaded before training by Python, even though the in-memory optional stage runs later.
            let edges = result
                .diagnostics
                .iter()
                .find(|s| s.starts_with("Number of network edges:"));
            let mut edge_emitted = false;
            for s in &result.diagnostics {
                if s.starts_with("Number of network edges:") || s.starts_with("Total ") {
                    continue;
                }
                if !edge_emitted
                    && (s == "# Caution: Resampling is disabled"
                        || s == "Iter TrainEss TrainNon TestSet")
                {
                    edge_emitted = true;
                    if let Some(edge) = edges {
                        println!("{edge}");
                    }
                }
                println!("{s}");
            }
            if args.test && cfg.network.is_some() && !cfg.guide_level {
                output::write_netscore(writer(&suffix(&args.output, ".netscore"))?, &result)?;
            }
            dump(
                &args.dump_json,
                &serde_json::json!({"config":cfg,"result":result,"inputs":{"matrix":identity(&args.input)?,"essential":identity(&args.essential)?,"nonessential":identity(&args.nonessential)?},"provenance":provenance(),"network_status":network_status.as_ref().err().map(ToString::to_string)}),
            )?;
            network_status?;
            output::write_bf(writer(&args.output)?, &result, &cfg)?;
            let metadata = serde_json::json!({"config":cfg,"inputs":{"matrix":identity(&args.input)?,"essential":identity(&args.essential)?,"nonessential":identity(&args.nonessential)?,"network":args.network.as_ref().map(|p|identity(p)).transpose()?,"alignment":args.alignment.as_ref().map(|p|identity(p)).transpose()?},"provenance":provenance(),"ranking_column":output::default_pr_column(&cfg),"legacy_defects_preserved":cfg.network_profile==NetworkProfile::Legacy});
            serde_json::to_writer_pretty(
                writer(&suffix(&args.output, ".metadata.json"))?,
                &metadata,
            )?;
        }
        Command::Pr(args) => {
            let ess = input::read_references(File::open(&args.essential)?)?;
            let non = input::read_references(File::open(&args.nonessential)?)?;
            let mut reader = csv_reader(&args.input)?;
            let header = reader.headers()?.clone();
            let column = args.column.unwrap_or_else(|| "BF".into());
            let index = header.iter().position(|x| x == column).ok_or_else(|| {
                Error::Invalid("Error! the column name is not in the file".into())
            })?;
            let mut genes = vec![];
            let mut scores = vec![];
            for row in reader.records() {
                let row = row?;
                genes.push(row[0].to_string());
                let value = &row[index];
                scores.push(
                    if value.is_empty() || value == "NA" || value == "NaN" || value == "nan" {
                        f64::NAN
                    } else {
                        value
                            .parse::<f64>()
                            .map_err(|_| Error::Invalid("invalid score".into()))?
                    },
                );
            }
            // Python opens/writes the header before a missing-essential failure.
            let mut w = writer(&args.output)?;
            match run_pr(&genes, &scores, &ess, &non) {
                Ok(rows) => {
                    output::write_pr(&mut w, &rows, &column)?;
                    dump(
                        &args.dump_json,
                        &serde_json::json!({"rows":rows,"inputs":{"scores":identity(&args.input)?,"essential":identity(&args.essential)?,"nonessential":identity(&args.nonessential)?},"provenance":provenance()}),
                    )?;
                    serde_json::to_writer_pretty(
                        writer(&suffix(&args.output, ".metadata.json"))?,
                        &serde_json::json!({"workflow":"pr","score_column":column,"inputs":{"scores":identity(&args.input)?,"essential":identity(&args.essential)?,"nonessential":identity(&args.nonessential)?},"provenance":provenance()}),
                    )?;
                }
                Err(err) => {
                    use std::io::Write;
                    writeln!(w, "Gene\t{column}\tRecall\tPrecision\tFDR")?;
                    return Err(err);
                }
            }
        }
    }
    Ok(())
}
fn csv_reader(p: &Path) -> Result<csv::Reader<File>> {
    Ok(csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(File::open(p)?))
}
fn main() {
    let _profile_session = core::profiling::Session::new();
    let _profile_process = core::profiling::span("cli_process");
    let args = std::env::args_os()
        .map(|a| match a.to_str() {
            Some("-Np") => "--pseudo-count".into(),
            Some("-NS") => "--no-resampling".into(),
            Some("-NB") => "--bootstrap-iterations".into(),
            Some("-m0") => "--loci-without-mismatch".into(),
            Some("-m1") => "--loci-with-mismatch".into(),
            _ => a,
        })
        .collect::<Vec<_>>();
    let cli = Cli::parse_from(args);
    if let Err(error) = run(cli.command) {
        eprintln!("{error}");
        std::process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bf_cli_defaults_to_numpy_and_legacy_with_explicit_native_opt_in() {
        let base = [
            "rsbagel", "bf", "-i", "fc", "-o", "bf", "-e", "ess", "-n", "non", "-c", "1",
        ];
        let Command::Bf(args) = Cli::try_parse_from(base).unwrap().command else {
            panic!("expected BF arguments");
        };
        assert_eq!(args.rng, "numpy");
        assert_eq!(args.profile, "legacy");
        assert_eq!(args.threads, core::workers::auto_workers());
        let native = base.into_iter().chain(["--rng-compat", "rust"]);
        let Command::Bf(args) = Cli::try_parse_from(native).unwrap().command else {
            panic!("expected BF arguments");
        };
        assert_eq!(args.rng, "rust");
    }

    #[test]
    fn workers_override_alias_and_provenance() {
        let base = [
            "rsbagel", "bf", "-i", "fc", "-o", "bf", "-e", "ess", "-n", "non", "-c", "1",
        ];
        for flag in ["--workers", "--threads"] {
            let Command::Bf(args) = Cli::try_parse_from(base.into_iter().chain([flag, "1"]))
                .unwrap()
                .command
            else {
                panic!("expected BF");
            };
            assert_eq!(args.threads, 1);
            let cfg = BfConfig {
                threads: args.threads,
                ..Default::default()
            };
            assert_eq!(serde_json::to_value(cfg).unwrap()["threads"], 1);
        }
        for value in ["0", "-1", "bad"] {
            assert!(Cli::try_parse_from(base.into_iter().chain(["--workers", value])).is_err());
        }
    }
}
