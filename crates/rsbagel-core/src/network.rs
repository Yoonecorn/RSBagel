//! Original BAGEL2 network model plus the single bounded opt-in correction.
//! Copyright (c) 2020 Hart Lab; MIT, see docs/BAGEL2-LICENSE.txt.
use crate::bf_types::{
    BfConfig, BfResult, NetworkEvent, NetworkProfile, NetworkTrace, Partition, Resampling,
};
use crate::kde::{arange, linregress, mean, round_hundredth, std, Kde};
use crate::resampling::{build_plans, validate_partition, Training};
use crate::{Error, Result};
use std::collections::{HashMap, HashSet};

pub fn fibonacci_weighted_sum(values: &[f64]) -> f64 {
    let mut value = 0.0;
    let mut p1 = 0.0;
    let mut c = 1.0;
    for &v in values {
        value += v / c;
        let p2 = p1;
        p1 = c;
        c = p1 + p2;
    }
    value
}

fn fit_network(
    result: &BfResult,
    partition: &Partition,
    essentials: &HashSet<&str>,
    nonessentials: &HashSet<&str>,
    check_cancel: &dyn Fn() -> Result<()>,
) -> Result<NetworkTrace> {
    check_cancel()?;
    let selected = |reference: &HashSet<&str>| -> Vec<f64> {
        partition
            .train
            .iter()
            .filter(|&&i| reference.contains(result.genes[i].gene.as_str()))
            .filter_map(|&i| result.network_scores[i])
            .collect()
    };
    let essential = selected(essentials);
    let nonessential = selected(nonessentials);
    let kess = Kde::new(&essential)?;
    let knon = Kde::new(&nonessential)?;
    let essential_mean = mean(&essential);
    let nonessential_mean = mean(&nonessential);
    let min = |x: &[f64]| x.iter().copied().fold(f64::INFINITY, f64::min);
    let max = |x: &[f64]| x.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let mut xmin = f64::INFINITY;
    let mut xmax = f64::NEG_INFINITY;
    for (step, score) in arange(max(&essential), min(&essential), -0.01)?
        .into_iter()
        .enumerate()
    {
        if step % 256 == 0 {
            check_cancel()?;
        }
        let e = kess.evaluate(score);
        let n = knon.evaluate(score);
        if e == 0.0 || n == 0.0 {
            continue;
        }
        if (e / n).log2() > -5.0 && score < essential_mean {
            xmin = xmin.min(score);
        }
    }
    for (step, score) in arange(min(&nonessential), max(&nonessential), 0.01)?
        .into_iter()
        .enumerate()
    {
        if step % 256 == 0 {
            check_cancel()?;
        }
        let e = kess.evaluate(score);
        let n = knon.evaluate(score);
        if e == 0.0 || n == 0.0 {
            continue;
        }
        if (e / n).log2() < 5.0 && score > nonessential_mean {
            xmax = xmax.max(score);
        }
    }
    let mut x = Vec::new();
    let mut y = Vec::new();
    for &i in &partition.train {
        if let Some(score) = result.network_scores[i] {
            if score >= xmin && score <= xmax {
                x.push(round_hundredth(score));
                y.push((kess.evaluate(score) / knon.evaluate(score)).log2());
            }
        }
    }
    let (slope, intercept) = linregress(&x, &y)?;
    Ok(NetworkTrace {
        partition: partition.clone(),
        essential,
        nonessential,
        xmin,
        xmax,
        regression_x: x,
        regression_y: y,
        slope,
        intercept,
    })
}

/// The caller retains `result` on an error, exposing precisely the completed
/// network partitions/events required for legacy partial-output diagnostics.
pub fn apply_network(
    result: &mut BfResult,
    essentials: &[String],
    nonessentials: &[String],
    config: &BfConfig,
) -> Result<()> {
    apply_network_with_progress(result, essentials, nonessentials, config, None)
}
pub fn apply_network_with_progress(
    result: &mut BfResult,
    essentials: &[String],
    nonessentials: &[String],
    config: &BfConfig,
    progress: Option<&(dyn Fn(usize, usize) -> std::ops::ControlFlow<()> + Sync)>,
) -> Result<()> {
    let Some(edges) = &config.network else {
        return Ok(());
    };
    if config.guide_level {
        return Ok(());
    }
    let _profile_span = crate::profiling::span(match config.network_profile {
        NetworkProfile::Legacy => "network_legacy",
        NetworkProfile::Corrected => "network_corrected",
    });
    let n = result.genes.len();
    let names: Vec<_> = result.genes.iter().map(|g| g.gene.clone()).collect();
    let indices: HashMap<&str, usize> = names
        .iter()
        .enumerate()
        .map(|(i, g)| (g.as_str(), i))
        .collect();
    let mut neighbors = vec![Vec::new(); n];
    let mut unique = vec![HashSet::new(); n];
    let mut edge_count = 0;
    for (a, b) in edges {
        if let (Some(&ia), Some(&ib)) = (indices.get(a.as_str()), indices.get(b.as_str())) {
            edge_count += 1;
            if unique[ia].insert(ib) {
                neighbors[ia].push(ib);
            }
            if unique[ib].insert(ia) {
                neighbors[ib].push(ia);
            }
        }
    }
    result
        .diagnostics
        .push(format!("Number of network edges: {edge_count}"));
    result
        .diagnostics
        .push("\nNetwork score calculation start\n".into());
    result.network_scores = neighbors
        .iter()
        .map(|adjacent| {
            if adjacent.is_empty() {
                None
            } else {
                let mut scores: Vec<_> =
                    adjacent.iter().map(|&i| result.genes[i].raw_mean).collect();
                // Any NaN neighbor makes the upstream weighted sum NaN regardless
                // of its sort position. Avoid an inconsistent sort comparator;
                // the resulting nonfinite training score still fails KDE.
                if scores.iter().any(|x| x.is_nan()) {
                    return Some(f64::NAN);
                }
                scores.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
                Some(fibonacci_weighted_sum(&scores))
            }
        })
        .collect();
    result.network_partitions.clear();
    result.network_events.clear();
    result.network_traces.clear();
    let ess: HashSet<_> = essentials.iter().map(String::as_str).collect();
    let non: HashSet<_> = nonessentials.iter().map(String::as_str).collect();
    let corrected = config.network_profile == NetworkProfile::Corrected;
    let count = result.partitions.len();
    let all: Vec<_> = (0..n).collect();
    let mut training = Training::new(&all, 10)?;
    let mut sampler = if corrected || config.replay.is_some() {
        None
    } else {
        Some(build_plans(&names, essentials, nonessentials, config)?.1)
    };
    let mut raw_boosted = vec![Vec::new(); n];
    let mut norm_boosted = vec![Vec::new(); n];
    let mut applied = vec![0usize; n];
    let mut contributions = vec![0usize; n];
    for iteration in 0..count {
        let check_cancel = || {
            if progress.is_some_and(|p| p(iteration, count).is_break()) {
                Err(Error::Cancelled)
            } else {
                Ok(())
            }
        };
        check_cancel()?;
        let partition = if corrected {
            result.partitions[iteration].clone()
        } else if let Some(replay) = &config.replay {
            match replay.network_iterations.get(iteration) {
                Some(p) => {
                    validate_partition(p, n, iteration, false)?;
                    p.clone()
                }
                None if config.resampling == Resampling::CrossValidation && iteration >= 10 => {
                    return Err(Error::Legacy("ValueError: high <= 0".into()))
                }
                None => {
                    return Err(Error::Invalid(
                        "legacy network replay requires separate network_iterations".into(),
                    ))
                }
            }
        } else {
            training.next_partition(
                sampler.as_mut().expect("sampler exists outside replay"),
                config.resampling != Resampling::CrossValidation,
                iteration,
            )?
        };
        result.network_partitions.push(partition.clone());
        let trace = fit_network(result, &partition, &ess, &non, &check_cancel).map_err(|e| {
            if corrected && !matches!(e, Error::Cancelled) {
                Error::Numerical(format!("corrected network model invalid: {e}"))
            } else {
                e
            }
        })?;
        if corrected && (!trace.slope.is_finite() || !trace.intercept.is_finite()) {
            return Err(Error::Numerical(
                "corrected network model invalid: non-finite regression".into(),
            ));
        }
        for &i in &partition.test {
            let score = result.network_scores[i];
            let predicted = score.map(|x| trace.slope * x + trace.intercept);
            let usable = predicted.is_some_and(f64::is_finite);
            let contribution = if corrected {
                predicted.filter(|x| x.is_finite()).unwrap_or(0.0)
            } else {
                predicted.unwrap_or(0.0)
            };
            contributions[i] += 1;
            if usable {
                applied[i] += 1;
            }
            raw_boosted[i].push(result.genes[i].raw_mean + contribution);
            if let Some(normalized) = result.genes[i].normalized {
                if corrected {
                    norm_boosted[i].push(normalized + contribution);
                } else {
                    raw_boosted[i].push(normalized + contribution);
                }
            }
            if config.dump_stages {
                if let Some(score) = score {
                    result.network_events.push(NetworkEvent {
                        gene: result.genes[i].gene.clone(),
                        score,
                        contribution,
                    });
                }
            }
        }
        if config.dump_stages {
            result.network_traces.push(trace);
        }
    }
    for (i, gene) in result.genes.iter_mut().enumerate() {
        if corrected && raw_boosted[i].is_empty() {
            // A valid model was fitted, but this individual gene received no
            // usable test contribution. This does not repair failed models.
            raw_boosted[i].push(gene.raw_mean);
            if let Some(normalized) = gene.normalized {
                norm_boosted[i].push(normalized);
            }
        }
        gene.boosted = Some(mean(&raw_boosted[i]));
        gene.boosted_std = Some(std(&raw_boosted[i]));
        if corrected && gene.normalized.is_some() {
            gene.boosted_normalized = Some(mean(&norm_boosted[i]));
        }
        gene.network_evidence = Some(applied[i] > 0);
        gene.network_contribution_count = Some(applied[i]);
        gene.network_zero_applied =
            Some(corrected && (applied[i] < contributions[i] || contributions[i] == 0));
    }
    if progress.is_some_and(|p| p(count, count).is_break()) {
        return Err(Error::Cancelled);
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn fibonacci_divisors() {
        assert_eq!(fibonacci_weighted_sum(&[1., 2., 4., 6., 10.]), 9.);
        assert_eq!(fibonacci_weighted_sum(&[]), 0.);
    }
    fn fixture(
        name: &str,
        no_resampling: bool,
        folds: usize,
    ) -> (
        BfResult,
        BfConfig,
        Vec<String>,
        Vec<String>,
        serde_json::Value,
    ) {
        use crate::bf_types::{GeneScore, RngMode};
        use std::path::Path;
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/network");
        let records: serde_json::Value = serde_json::from_slice(
            &std::fs::read(root.join(format!("{name}.trace.json"))).unwrap(),
        )
        .unwrap();
        let last = records.as_array().unwrap().last().unwrap().clone();
        let raw = last["bf_mean"].as_object().unwrap();
        let genes: Vec<_> = raw
            .iter()
            .map(|(gene, value)| GeneScore {
                gene: gene.clone(),
                raw_mean: value.as_f64().unwrap(),
                mean: value.as_f64().unwrap(),
                std: 0.,
                count: 1,
                normalized: last["bf_norm"][gene].as_f64(),
                boosted: None,
                boosted_std: None,
                boosted_normalized: None,
                network_evidence: None,
                network_contribution_count: None,
                network_zero_applied: None,
            })
            .collect();
        let names: Vec<_> = genes.iter().map(|g| g.gene.clone()).collect();
        let read_refs = |name: &str| {
            std::fs::read_to_string(root.join(name))
                .unwrap()
                .lines()
                .skip(1)
                .map(|v| v.split('\t').next().unwrap().to_string())
                .collect::<Vec<_>>()
        };
        let ess = read_refs("essential.tsv");
        let non = read_refs("nonessential.tsv");
        let edges = std::fs::read_to_string(root.join("network.tsv"))
            .unwrap()
            .lines()
            .map(|l| {
                let mut f = l.split('\t');
                (f.next().unwrap().to_string(), f.next().unwrap().to_string())
            })
            .collect();
        let config = BfConfig {
            columns: vec![0, 1],
            no_resampling,
            folds,
            seed: 5411,
            rng_mode: RngMode::Numpy,
            network: Some(edges),
            dump_stages: true,
            ..BfConfig::default()
        };
        let partitions = build_plans(&names, &ess, &non, &config).unwrap().0;
        let result = BfResult {
            genes,
            guides: Vec::new(),
            gene_guides: Vec::new(),
            partitions,
            guide_iterations: Vec::new(),
            traces: Vec::new(),
            network_partitions: Vec::new(),
            network_scores: Vec::new(),
            network_events: Vec::new(),
            network_traces: Vec::new(),
            diagnostics: Vec::new(),
            correction_coefficients: None,
            correction_trace: None,
        };
        (result, config, ess, non, last)
    }
    #[test]
    fn pinned_legacy_ns_network_nonfinite_and_values() {
        let (mut result, config, ess, non, reference) = fixture("ns_network", true, 10);
        apply_network(&mut result, &ess, &non, &config).unwrap();
        assert_eq!(result.network_partitions.len(), 1);
        assert_eq!(
            result
                .genes
                .iter()
                .filter(|g| g.boosted.unwrap().is_nan())
                .count(),
            108
        );
        for gene in &result.genes {
            let observed = reference["boostedbf"][&gene.gene].as_array().unwrap();
            if observed.is_empty() {
                assert!(gene.boosted.unwrap().is_nan());
            } else {
                assert!(
                    (gene.boosted.unwrap() - observed[0].as_f64().unwrap()).abs() < 1e-8,
                    "{}",
                    gene.gene
                );
            }
        }
    }
    #[test]
    fn pinned_legacy_equalization_double_append() {
        let (mut result, mut config, ess, non, reference) = fixture("equal4_network", false, 10);
        config.equalise_guides = Some(4.);
        apply_network(&mut result, &ess, &non, &config).unwrap();
        for gene in &result.genes {
            let values = reference["boostedbf"][&gene.gene].as_array().unwrap();
            let expected =
                values.iter().map(|v| v.as_f64().unwrap()).sum::<f64>() / values.len() as f64;
            assert!(
                (gene.boosted.unwrap() - expected).abs() < 1e-8,
                "{} {} {}",
                gene.gene,
                gene.boosted.unwrap(),
                expected
            );
            assert!(gene.boosted_normalized.is_none());
        }
    }
    #[test]
    fn corrected_reuses_partition_and_keeps_separate_values() {
        let (mut result, mut config, ess, non, _) = fixture("equal4_network", false, 11);
        config.network_profile = NetworkProfile::Corrected;
        config.equalise_guides = Some(4.);
        apply_network(&mut result, &ess, &non, &config).unwrap();
        assert_eq!(result.network_partitions.len(), 11);
        for (main, network) in result.partitions.iter().zip(&result.network_partitions) {
            assert_eq!(main.train, network.train);
            assert_eq!(main.test, network.test);
        }
        for gene in &result.genes {
            let raw_addition = gene.boosted.unwrap() - gene.raw_mean;
            let normalized_addition = gene.boosted_normalized.unwrap() - gene.normalized.unwrap();
            assert!((raw_addition - normalized_addition).abs() < 1e-12);
            assert!(gene.network_evidence.unwrap());
        }
    }
    #[test]
    fn corrected_missing_gene_evidence_is_zero_but_invalid_model_fails() {
        let (mut result, mut config, ess, non, _) = fixture("ns_network", true, 10);
        config.network_profile = NetworkProfile::Corrected;
        let missing = "G119";
        config
            .network
            .as_mut()
            .unwrap()
            .retain(|(a, b)| a != missing && b != missing);
        apply_network(&mut result, &ess, &non, &config).unwrap();
        assert!(result.genes.iter().all(|g| g.boosted.unwrap().is_finite()));
        let gene = result.genes.iter().find(|g| g.gene == missing).unwrap();
        assert_eq!(gene.boosted, Some(gene.raw_mean));
        assert_eq!(gene.network_evidence, Some(false));
        assert_eq!(gene.network_zero_applied, Some(true));
        let (mut result, mut config, ess, non, _) = fixture("ns_network", true, 10);
        config.network_profile = NetworkProfile::Corrected;
        config.network = Some(vec![]);
        assert!(apply_network(&mut result, &ess, &non, &config)
            .unwrap_err()
            .to_string()
            .contains("corrected network model invalid"));
        assert!(result.genes.iter().all(|g| g.boosted.is_none()));
    }
    #[test]
    fn legacy_nondefault_cv_retains_completed_partial_trace() {
        let (mut result, config, ess, non, _) = fixture("cv_network", false, 11);
        let err = apply_network(&mut result, &ess, &non, &config).unwrap_err();
        assert_eq!(err.to_string(), "ValueError: high <= 0");
        assert_eq!(result.network_partitions.len(), 10);
        assert_eq!(result.network_events.len(), 120);
    }
    #[test]
    fn cancellation_inside_network_grid_preserves_raw_scores() {
        use std::sync::atomic::{AtomicUsize, Ordering};
        for profile in [NetworkProfile::Legacy, NetworkProfile::Corrected] {
            let (mut result, mut cfg, ess, non, _) = fixture("ns_network", true, 10);
            cfg.network_profile = profile;
            let raw: Vec<_> = result.genes.iter().map(|g| g.raw_mean.to_bits()).collect();
            let calls = AtomicUsize::new(0);
            let callback = |_, _| {
                if calls.fetch_add(1, Ordering::Relaxed) >= 2 {
                    std::ops::ControlFlow::Break(())
                } else {
                    std::ops::ControlFlow::Continue(())
                }
            };
            assert!(matches!(
                apply_network_with_progress(&mut result, &ess, &non, &cfg, Some(&callback)),
                Err(Error::Cancelled)
            ));
            assert_eq!(
                raw,
                result
                    .genes
                    .iter()
                    .map(|g| g.raw_mean.to_bits())
                    .collect::<Vec<_>>()
            );
            assert!(result.network_events.is_empty());
        }
    }
}
