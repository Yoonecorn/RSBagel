//! Exploratory pure-memory pipeline, with explicit source precision bridges.
use rsbagel_core as core;
use std::fs::File;
use std::io::BufWriter;
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let _session = core::profiling::Session::new();
    let _pipeline = core::profiling::span("pure_memory_pipeline");
    let args: Vec<_> = std::env::args().collect();
    if args.len() != 6 {
        return Err(
            "usage: profile_memory COUNTS ESSENTIAL NONESSENTIAL BF_CONFIG_JSON OUTPUT_DIRECTORY"
                .into(),
        );
    }
    let counts = core::input::read_matrix(File::open(&args[1])?, true)?;
    let ess = core::input::read_references(File::open(&args[2])?)?;
    let non = core::input::read_references(File::open(&args[3])?)?;
    let mut cfg: core::BfConfig = serde_json::from_reader(File::open(&args[4])?)?;
    if cfg.rng_mode != core::RngMode::Numpy || cfg.replay.is_some() {
        return Err("profile_memory requires production NumPy sampling without replay".into());
    }
    let fc_cfg = core::FcConfig {
        controls: vec![0],
        control_labels: vec![counts.columns[0].clone()],
        ..Default::default()
    };
    let fc = core::run_fc(&counts, &fc_cfg)?;
    let bf_input = core::fc_to_bf_input(&fc.fold_change)?;
    cfg.column_labels = cfg
        .columns
        .iter()
        .map(|&c| bf_input.columns[c].clone())
        .collect();
    let bf = core::run_bf(&bf_input, &ess, &non, &cfg, None)?;
    let column = core::output::default_pr_column(&cfg);
    let (genes, scores) = core::bf_to_pr_input(&bf, &cfg, column)?;
    let pr = core::run_pr(&genes, &scores, &ess, &non)?;
    let directory = std::path::Path::new(&args[5]);
    core::output::write_matrix(
        BufWriter::new(File::create(directory.join("result.foldchange"))?),
        &fc.fold_change,
        false,
    )?;
    core::output::write_matrix(
        BufWriter::new(File::create(directory.join("result.normed_readcount"))?),
        &fc.normalized,
        true,
    )?;
    core::output::write_bf(
        BufWriter::new(File::create(directory.join("result.bf"))?),
        &bf,
        &cfg,
    )?;
    core::output::write_pr(
        BufWriter::new(File::create(directory.join("result.pr"))?),
        &pr,
        column,
    )?;
    Ok(())
}
