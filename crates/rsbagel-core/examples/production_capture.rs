//! Compact default-profile HAP1 bootstrap1000 qualification capture.
//! This diagnostic executable is not the production CLI or a benchmark harness.
use rsbagel_core::{input, output, run_bf, BfConfig, Resampling, SOURCE_COMMIT, SOURCE_DIRTY};
use serde_json::json;
use sha2::{Digest, Sha256};
use std::{
    fs::{self, File, OpenOptions},
    io::BufWriter,
    path::PathBuf,
};

fn create(path: PathBuf) -> std::io::Result<BufWriter<File>> {
    Ok(BufWriter::new(
        OpenOptions::new().write(true).create_new(true).open(path)?,
    ))
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let arguments: Vec<_> = std::env::args_os().skip(1).collect();
    if arguments.len() != 5 {
        return Err("usage: production_capture FC ESS NON OUTPUT_DIRECTORY WORKERS".into());
    }
    let matrix = input::read_matrix(File::open(&arguments[0])?, false)?;
    let essential = input::read_references(File::open(&arguments[1])?)?;
    let nonessential = input::read_references(File::open(&arguments[2])?)?;
    let directory = PathBuf::from(&arguments[3]);
    fs::create_dir(&directory)?;
    let config = BfConfig {
        columns: (0..matrix.columns.len()).collect(),
        column_labels: matrix.columns.clone(),
        resampling: Resampling::Bootstrap,
        bootstrap_iterations: 1000,
        threads: arguments[4].to_str().ok_or("invalid workers")?.parse()?,
        // RNG, seed and network use the production defaults under qualification.
        ..BfConfig::default()
    };
    serde_json::to_writer_pretty(
        create(directory.join("candidate.json"))?,
        &json!({
            "source_commit": SOURCE_COMMIT, "source_dirty": SOURCE_DIRTY, "config": config,
            "precision": "f64", "scope": "default-profile pure-memory bootstrap1000 capture",
            "native_rng_runs": 0
        }),
    )?;
    let result = run_bf(&matrix, &essential, &nonessential, &config, None)?;
    let genes: Vec<_> = result.genes.iter().map(|gene| &gene.gene).collect();
    let values: Vec<_> = result
        .genes
        .iter()
        .map(|gene| [gene.raw_mean, gene.std, gene.count as f64])
        .collect();
    let partitions: Vec<_> = result
        .partitions
        .iter()
        .map(|partition| {
            let mut hash = Sha256::new();
            for value in [
                partition.iteration,
                partition.train.len(),
                partition.test.len(),
            ]
            .into_iter()
            .chain(partition.train.iter().copied())
            .chain(partition.test.iter().copied())
            {
                hash.update((value as u64).to_le_bytes());
            }
            json!({"iteration": partition.iteration, "train_count": partition.train.len(),
            "test_count": partition.test.len(), "sha256": format!("{:x}", hash.finalize())})
        })
        .collect();
    serde_json::to_writer(create(directory.join("genes.json"))?, &genes)?;
    serde_json::to_writer(create(directory.join("mean-std-count.json"))?, &values)?;
    serde_json::to_writer_pretty(
        create(directory.join("partition-hashes.json"))?,
        &partitions,
    )?;
    output::write_bf(create(directory.join("result.bf"))?, &result, &config)?;
    Ok(())
}
