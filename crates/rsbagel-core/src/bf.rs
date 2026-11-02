//! Faithful main-stage BAGEL2 fit and scoring; optional corrections are separate.
use crate::bf_types::*;
use crate::kde::{arange, linregress, mean, round_hundredth, std, Kde};
use crate::resampling::build_plans;
use crate::{Error, Matrix, Result};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::ops::ControlFlow;

type Progress<'a> = Option<&'a (dyn Fn(usize, usize) -> ControlFlow<()> + Sync)>;
struct Iteration {
    trace: Option<IterationTrace>,
    scores: Vec<(usize, Vec<f64>)>,
}

pub fn run_bf_base(
    model: &Matrix,
    essentials: &[String],
    nonessentials: &[String],
    config: &BfConfig,
    progress: Progress<'_>,
) -> Result<BfResult> {
    let _profile_span = crate::profiling::span("bf_total");
    model.validate()?;
    if [config.equalise_guides, config.equalise_reps]
        .into_iter()
        .flatten()
        .any(|v| !v.is_finite() || v.fract() != 0.0)
    {
        return Err(Error::Invalid(
            "equalization counts must be finite integers".into(),
        ));
    }
    if config.columns.is_empty() || config.columns.iter().any(|&c| c >= model.cols()) {
        return Err(Error::Invalid("Invalid columns".into()));
    }
    if config.threads == 0 {
        return Err(Error::Invalid(
            "worker count must be greater than zero".into(),
        ));
    }
    if let Some(records) = &config.alignment {
        let mut alignment_ids = HashSet::new();
        if records
            .iter()
            .any(|a| !alignment_ids.insert(a.guide.as_str()))
        {
            return Err(Error::Legacy("Please check align-info file".into()));
        }
    }
    let alignment: HashMap<&str, &AlignmentRecord> = config
        .alignment
        .as_deref()
        .unwrap_or_default()
        .iter()
        .map(|a| (a.guide.as_str(), a))
        .collect();
    let rows: Vec<usize> = (0..model.rows())
        .filter(|&r| {
            alignment.get(model.row_ids[r].as_str()).is_none_or(|a| {
                a.perfect_loci <= config.max_perfect_loci
                    && a.mismatch_loci <= config.max_mismatch_loci
            })
        })
        .collect();
    let discarded = model.rows() - rows.len();
    let mut guide_ids = HashSet::new();
    let mut genes = Vec::new();
    let mut gene_indices = HashMap::new();
    let mut gene_guides: Vec<Vec<usize>> = Vec::new();
    let mut guide_gene = Vec::new();
    for (guide, &row) in rows.iter().enumerate() {
        if !guide_ids.insert(model.row_ids[row].as_str()) {
            return Err(Error::Legacy("Error! sgRNA tag duplicates".into()));
        }
        let gene = model.row_genes[row].as_str();
        let index = if let Some(&i) = gene_indices.get(gene) {
            i
        } else {
            let i = genes.len();
            genes.push(gene.to_owned());
            gene_indices.insert(gene, i);
            gene_guides.push(Vec::new());
            i
        };
        gene_guides[index].push(guide);
        guide_gene.push(index);
    }
    let essential: HashSet<&str> = essentials.iter().map(String::as_str).collect();
    let nonessential: HashSet<&str> = nonessentials.iter().map(String::as_str).collect();
    let threshold = 2.0_f64.powf(
        -1.1535
            * ((genes
                .iter()
                .filter(|g| nonessential.contains(g.as_str()))
                .count() as f64)
                + 13.324)
                .ln()
            + 0.7728,
    );
    let mut unique_columns = Vec::new();
    let mut seen_columns = HashSet::new();
    let mut unique_positions = Vec::new();
    for (position, &c) in config.columns.iter().enumerate() {
        if seen_columns.insert(c) {
            unique_columns.push(c);
            unique_positions.push(position);
        }
    }
    let (partitions, _sampler) = build_plans(&genes, essentials, nonessentials, config)?;
    let total = partitions.len();
    if progress.is_some_and(|f| f(0, total).is_break()) {
        return Err(Error::Cancelled);
    }
    let fit_one = |partition: &Partition| -> Result<Iteration> {
        iteration(
            model,
            &rows,
            &genes,
            &gene_guides,
            &essential,
            &nonessential,
            &unique_columns,
            &config.columns,
            partition,
            threshold,
            config.dump_stages,
        )
    };
    // A fixed indexed collect and subsequent serial reduction preserves output order.
    // Dispatch is not an optimization experiment; no profile-dependent math is used.
    let results: Vec<Result<Iteration>> = if config.threads == 1 {
        let mut result = Vec::with_capacity(total);
        for (i, p) in partitions.iter().enumerate() {
            result.push(fit_one(p));
            if result.last().is_some_and(Result::is_err) {
                break;
            }
            if progress.is_some_and(|f| f(i + 1, total).is_break()) {
                return Err(Error::Cancelled);
            }
        }
        result
    } else {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(config.threads)
            .build()
            .map_err(|e| Error::Invalid(e.to_string()))?;
        let mut result = Vec::with_capacity(total);
        for chunk in partitions.chunks(config.threads) {
            let batch: Vec<Result<Iteration>> =
                pool.install(|| chunk.par_iter().map(fit_one).collect());
            for item in batch {
                let failed = item.is_err();
                result.push(item);
                if failed {
                    break;
                }
                if progress.is_some_and(|f| f(result.len(), total).is_break()) {
                    return Err(Error::Cancelled);
                }
            }
            if result.last().is_some_and(Result::is_err) {
                break;
            }
        }
        result
    };
    let _profile_aggregation = crate::profiling::span("guide_gene_aggregation");
    let mut guide_iterations = vec![Vec::new(); rows.len()];
    let mut traces = Vec::new();
    for result in results {
        let iter = result?;
        for (guide, scores) in iter.scores {
            guide_iterations[guide].push(scores);
        }
        if let Some(trace) = iter.trace {
            traces.push(trace);
        }
    }
    let mut guides = Vec::with_capacity(rows.len());
    for (guide, &row) in rows.iter().enumerate() {
        let observations = &guide_iterations[guide];
        if observations.is_empty() && (config.guide_level || config.alignment.is_some()) {
            return Err(Error::Legacy("IndexError: list index out of range".into()));
        }
        let replicate_values: Vec<Vec<f64>> = (0..config.columns.len())
            .map(|rep| observations.iter().map(|v| v[rep]).collect())
            .collect();
        let means: Vec<f64> = replicate_values.iter().map(|v| mean(v)).collect();
        let deviations = replicate_values.iter().map(|v| std(v)).collect();
        let total = unique_positions.iter().fold(0.0, |sum, &p| sum + means[p]);
        guides.push(GuideScore {
            guide: model.row_ids[row].clone(),
            gene: genes[guide_gene[guide]].clone(),
            mean_replicates: means,
            std_replicates: deviations,
            total,
            count: observations.len(),
        });
    }
    let gene_scores = genes
        .iter()
        .enumerate()
        .map(|(g, name)| {
            let count = guide_iterations[gene_guides[g][0]].len();
            let values: Vec<f64> = (0..count)
                .map(|i| {
                    gene_guides[g].iter().fold(0.0, |sum, &guide| {
                        sum + guide_iterations[guide][i].iter().fold(0.0, |r, v| r + v)
                    })
                })
                .collect();
            let raw_mean = mean(&values);
            let normalized = if !config.guide_level {
                config
                    .equalise_guides
                    .filter(|v| *v != 0.0)
                    .map(|n| raw_mean * (n / gene_guides[g].len() as f64))
            } else {
                None
            };
            GeneScore {
                gene: name.clone(),
                raw_mean,
                mean: raw_mean,
                std: std(&values),
                count,
                normalized,
                boosted: None,
                boosted_std: None,
                boosted_normalized: None,
                network_evidence: None,
                network_contribution_count: None,
                network_zero_applied: None,
            }
        })
        .collect();
    let mut diagnostics = Vec::new();
    if config.alignment.is_some() {
        diagnostics.push(format!(
            "Total {discarded} multi-targeting gRNAs are discarded"
        ));
    }
    diagnostics.push(format!("Number of unique genes:  {}", genes.len()));
    diagnostics.push(format!(
        "Number of reference essentials: {}",
        essentials.len()
    ));
    diagnostics.push(format!(
        "Number of reference nonessentials: {}",
        nonessentials.len()
    ));
    if config.no_resampling {
        diagnostics.push("# Caution: Resampling is disabled".into());
    }
    diagnostics.push("Iter TrainEss TrainNon TestSet".into());
    for p in &partitions {
        diagnostics.push(format!(
            "{} {} {} {}",
            p.iteration,
            p.train
                .iter()
                .filter(|&&g| essential.contains(genes[g].as_str()))
                .count(),
            p.train
                .iter()
                .filter(|&&g| nonessential.contains(genes[g].as_str()))
                .count(),
            p.test.len()
        ));
    }
    Ok(BfResult {
        genes: gene_scores,
        guides,
        gene_guides,
        partitions,
        guide_iterations,
        traces,
        network_partitions: Vec::new(),
        network_scores: Vec::new(),
        diagnostics,
        correction_coefficients: None,
        correction_trace: None,
        network_events: Vec::new(),
        network_traces: Vec::new(),
    })
}

#[allow(clippy::too_many_arguments)]
fn iteration(
    model: &Matrix,
    rows: &[usize],
    genes: &[String],
    gene_guides: &[Vec<usize>],
    essential: &HashSet<&str>,
    nonessential: &HashSet<&str>,
    training_columns: &[usize],
    scoring_columns: &[usize],
    partition: &Partition,
    threshold: f64,
    dump: bool,
) -> Result<Iteration> {
    let _profile_fit = crate::profiling::span("iteration_fit_and_kde");
    let observations = |reference: &HashSet<&str>| -> Vec<f64> {
        let mut values = Vec::new();
        for &gene in &partition.train {
            if reference.contains(genes[gene].as_str()) {
                for &guide in &gene_guides[gene] {
                    for &column in training_columns {
                        values.push(model.at(rows[guide], column));
                    }
                }
            }
        }
        values
    };
    let essential_values = observations(essential);
    let nonessential_values = observations(nonessential);
    let kess = Kde::new(&essential_values)?;
    let knon = Kde::new(&nonessential_values)?;
    let _profile_evaluation = crate::profiling::span("kde_evaluate_support_ratio_lookup");
    let grid = arange(-10.0, 2.0, 0.01)?;
    let density = knon.evaluate_many(&grid);
    let support: Vec<usize> = density
        .iter()
        .enumerate()
        .filter_map(|(i, &v)| (v > threshold).then_some(i))
        .collect();
    if support.is_empty() {
        return Err(Error::Legacy(
            "ValueError: min() arg is an empty sequence".into(),
        ));
    }
    let xmin = round_hundredth(grid[support[0]]);
    let subgrid = arange(xmin, grid[*support.last().unwrap()], 0.01)?;
    let logratios: Vec<f64> = subgrid
        .iter()
        .map(|&x| (kess.evaluate(x) / knon.evaluate(x)).log2())
        .collect();
    if logratios.is_empty() || logratios.iter().any(|v| v.is_nan()) {
        return Err(Error::Numerical(
            "invalid or empty KDE density-ratio fit grid".into(),
        ));
    }
    let min = logratios.iter().copied().fold(f64::INFINITY, f64::min);
    let minimum_indices: Vec<usize> = logratios
        .iter()
        .enumerate()
        .filter_map(|(i, &v)| (v == min).then_some(i))
        .collect();
    if minimum_indices.len() != 1 {
        return Err(Error::Legacy(
            "multiple exact density-ratio minima cannot be converted to a scalar bound".into(),
        ));
    }
    let xmax = round_hundredth(subgrid[minimum_indices[0]]);
    let lookup_grid = arange(xmin, xmax + 0.01, 0.01)?;
    let mut lookup = HashMap::new();
    for (index, x) in lookup_grid.into_iter().enumerate() {
        // These aranges normally share their start and effective step. Reuse only
        // an identical f64 evaluation point; a rounded key is not sufficient.
        let ratio = match subgrid.get(index) {
            Some(point) if point.to_bits() == x.to_bits() => logratios[index],
            _ => (kess.evaluate(x) / knon.evaluate(x)).log2(),
        };
        lookup.insert((x * 100.0).round_ties_even() as i64, ratio);
    }
    drop(_profile_evaluation);
    let mut regression_x = Vec::new();
    let mut regression_y = Vec::new();
    for &gene in &partition.train {
        for &guide in &gene_guides[gene] {
            for &column in training_columns {
                let fc = model.at(rows[guide], column);
                if fc >= xmin && fc <= xmax {
                    let key = (fc * 100.0).round_ties_even() as i64;
                    let y = *lookup
                        .get(&key)
                        .ok_or_else(|| Error::Legacy(format!("KeyError: {key}")))?;
                    regression_x.push(key as f64 / 100.0);
                    regression_y.push(y);
                }
            }
        }
    }
    let (slope, intercept) = linregress(&regression_x, &regression_y)
        .map_err(|_| Error::Legacy("Regression failed. Check quality of the screen".into()))?;
    drop(_profile_fit);
    let _profile_scoring = crate::profiling::span("iteration_bf_score");
    let mut scores = Vec::new();
    for &gene in &partition.test {
        for &guide in &gene_guides[gene] {
            scores.push((
                guide,
                scoring_columns
                    .iter()
                    .map(|&column| slope * model.at(rows[guide], column) + intercept)
                    .collect(),
            ));
        }
    }
    let trace = if dump {
        Some(IterationTrace {
            partition: partition.clone(),
            essential: essential_values,
            nonessential: nonessential_values,
            essential_factor: kess.factor,
            nonessential_factor: knon.factor,
            essential_covariance: kess.covariance,
            nonessential_covariance: knon.covariance,
            grid,
            nonessential_density: density,
            support_indices: support,
            subgrid,
            logratios,
            minimum_indices,
            xmin,
            xmax,
            regression_x,
            regression_y,
            slope,
            intercept,
        })
    } else {
        None
    };
    Ok(Iteration { trace, scores })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::{Path, PathBuf};
    fn root() -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
    }
    fn reference(path: &Path) -> Vec<String> {
        std::fs::read_to_string(path)
            .unwrap()
            .lines()
            .skip(1)
            .map(|l| l.split('\t').next().unwrap().to_owned())
            .collect()
    }
    fn load_matrix(path: &Path) -> Matrix {
        let text = std::fs::read_to_string(path).unwrap();
        let mut lines = text.lines();
        let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
        let columns = header.iter().skip(2).map(|v| (*v).to_owned()).collect();
        let mut m = Matrix {
            index_name: header[0].to_owned(),
            row_ids: Vec::new(),
            row_genes: Vec::new(),
            columns,
            values: Vec::new(),
        };
        for line in lines {
            let f: Vec<&str> = line.split('\t').collect();
            m.row_ids.push(f[0].to_owned());
            m.row_genes.push(f[1].to_owned());
            m.values
                .extend(f[2..].iter().map(|v| v.parse::<f64>().unwrap()));
        }
        m
    }
    fn npy(path: &Path) -> Vec<f64> {
        let bytes = std::fs::read(path).unwrap();
        assert_eq!(&bytes[..6], b"\x93NUMPY");
        let (offset, length) = if bytes[6] == 1 {
            (10, u16::from_le_bytes([bytes[8], bytes[9]]) as usize)
        } else {
            (
                12,
                u32::from_le_bytes(bytes[8..12].try_into().unwrap()) as usize,
            )
        };
        let header = std::str::from_utf8(&bytes[offset..offset + length]).unwrap();
        bytes[offset + length..]
            .as_chunks::<8>()
            .0
            .iter()
            .map(|b| {
                if header.contains("<f8") {
                    f64::from_le_bytes(*b)
                } else if header.contains("<i8") {
                    i64::from_le_bytes(*b) as f64
                } else {
                    panic!("unexpected NPY type: {header}")
                }
            })
            .collect()
    }
    fn close(name: &str, actual: &[f64], expected: &[f64], atol: f64, rtol: f64) {
        assert_eq!(actual.len(), expected.len(), "{name} shape");
        for (i, (&a, &b)) in actual.iter().zip(expected).enumerate() {
            assert!(
                (a - b).abs() <= atol + rtol * b.abs() || (a.is_nan() && b.is_nan()),
                "{name}[{i}]: {a:?} vs {b:?}, abs {}",
                (a - b).abs()
            );
        }
    }
    fn run_golden(name: &str, mode: Resampling) -> BfResult {
        let dir = if name.starts_with("hap1") {
            std::env::var_os("RSBAGEL_HAP1_GOLDEN_DIR")
                .map(PathBuf::from)
                .unwrap_or_else(|| root().join("tests/external/hap1-golden"))
        } else {
            root().join("tests/golden").join(name)
        };
        let m = load_matrix(&dir.join("reference.foldchange"));
        let (ess, non) = if name.starts_with("hap1") {
            (
                root().join("upstream/BAGEL2/CEGv2.txt"),
                root().join("upstream/BAGEL2/NEGv1.txt"),
            )
        } else {
            (
                dir.join("inputs/essentials.tsv"),
                dir.join("inputs/nonessentials.tsv"),
            )
        };
        let config = BfConfig {
            columns: (0..m.cols()).collect(),
            replay: Some(
                serde_json::from_str(&std::fs::read_to_string(dir.join("replay.json")).unwrap())
                    .unwrap(),
            ),
            dump_stages: true,
            rng_mode: RngMode::Numpy,
            resampling: mode,
            bootstrap_iterations: 100,
            ..Default::default()
        };
        let result = run_bf_base(&m, &reference(&ess), &reference(&non), &config, None).unwrap();
        for t in &result.traces {
            let p = format!("iter_{:04}_", t.partition.iteration);
            let a = dir.join("arrays");
            close(
                &format!("{p}essential"),
                &t.essential,
                &npy(&a.join(format!("{p}ess_train_fc_flat_list.npy"))),
                0.,
                0.,
            );
            close(
                &format!("{p}nonessential"),
                &t.nonessential,
                &npy(&a.join(format!("{p}non_train_fc_flat_list.npy"))),
                0.,
                0.,
            );
            close(
                &format!("{p}factor"),
                &[t.essential_factor, t.nonessential_factor],
                &[
                    npy(&a.join(format!("{p}kess_factor.npy")))[0],
                    npy(&a.join(format!("{p}knon_factor.npy")))[0],
                ],
                1e-10,
                1e-8,
            );
            close(
                &format!("{p}covariance"),
                &[t.essential_covariance, t.nonessential_covariance],
                &[
                    npy(&a.join(format!("{p}kess_covariance.npy")))[0],
                    npy(&a.join(format!("{p}knon_covariance.npy")))[0],
                ],
                1e-10,
                1e-8,
            );
            close(
                &format!("{p}grid"),
                &t.grid,
                &npy(&a.join(format!("{p}x.npy"))),
                0.,
                0.,
            );
            close(
                &format!("{p}density"),
                &t.nonessential_density,
                &npy(&a.join(format!("{p}nonfitx.npy"))),
                1e-10,
                1e-8,
            );
            close(
                &format!("{p}subgrid"),
                &t.subgrid,
                &npy(&a.join(format!("{p}subx.npy"))),
                0.,
                0.,
            );
            close(
                &format!("{p}support"),
                &t.support_indices
                    .iter()
                    .map(|i| *i as f64)
                    .collect::<Vec<_>>(),
                &npy(&a.join(format!("{p}density_support_indices.npy"))),
                0.,
                0.,
            );
            close(
                &format!("{p}minima"),
                &t.minimum_indices
                    .iter()
                    .map(|i| *i as f64)
                    .collect::<Vec<_>>(),
                &npy(&a.join(format!("{p}fit_minimum_indices.npy"))),
                0.,
                0.,
            );
            close(
                &format!("{p}regression_x"),
                &t.regression_x,
                &npy(&a.join(format!("{p}regression_x.npy"))),
                0.,
                0.,
            );
            close(
                &format!("{p}regression_y"),
                &t.regression_y,
                &npy(&a.join(format!("{p}regression_y.npy"))),
                1e-10,
                1e-8,
            );
            close(
                &format!("{p}coefficients"),
                &[t.slope, t.intercept],
                &npy(&a.join(format!("{p}regression_slope_intercept.npy"))),
                1e-9,
                1e-8,
            );
        }
        let a = dir.join("arrays");
        let scores: Vec<f64> = result
            .guide_iterations
            .iter()
            .flat_map(|v| v.iter().flat_map(|r| r.iter().copied()))
            .collect();
        close(
            "guide iteration scores",
            &scores,
            &npy(&a.join("guide_iteration_bf.npy")),
            1e-8,
            0.,
        );
        let summaries: Vec<f64> = result
            .genes
            .iter()
            .flat_map(|g| [g.raw_mean, g.std, g.count as f64])
            .collect();
        close(
            "gene summaries",
            &summaries,
            &npy(&a.join("gene_mean_std_count.npy")),
            1e-8,
            0.,
        );
        result
    }
    #[test]
    fn synthetic_replay_stages() {
        run_golden("synthetic-cv-final-a", Resampling::CrossValidation);
    }
    #[test]
    fn synthetic_bootstrap_replay_stages() {
        run_golden("synthetic-bootstrap-final", Resampling::Bootstrap);
    }
    #[test]
    fn synthetic_small_sample_replay_stages() {
        run_golden("synthetic-small-sample-final", Resampling::SmallSample);
    }
    #[test]
    #[ignore = "HAP1 regression fixture is not distributed in the snapshot"]
    fn hap1_replay_stages() {
        run_golden("hap1-cv-final", Resampling::CrossValidation);
    }
    #[test]
    fn generated_numpy_plans_and_thread_identity() {
        let dir = root().join("tests/golden/synthetic-cv-final-a");
        let m = load_matrix(&dir.join("reference.foldchange"));
        let ess = reference(&dir.join("inputs/essentials.tsv"));
        let non = reference(&dir.join("inputs/nonessentials.tsv"));
        let replay: ReplayManifest =
            serde_json::from_str(&std::fs::read_to_string(dir.join("replay.json")).unwrap())
                .unwrap();
        let mut c = BfConfig {
            columns: vec![0, 1],
            rng_mode: RngMode::Numpy,
            ..Default::default()
        };
        let result = run_bf_base(&m, &ess, &non, &c, None).unwrap();
        assert_eq!(
            serde_json::to_vec(&result.partitions).unwrap(),
            serde_json::to_vec(&replay.iterations).unwrap()
        );
        let serial = serde_json::to_vec(&result).unwrap();
        for threads in [2, 4, 8] {
            c.threads = threads;
            assert_eq!(
                serial,
                serde_json::to_vec(&run_bf_base(&m, &ess, &non, &c, None).unwrap()).unwrap(),
                "worker{threads}"
            );
        }
        for (name, mode) in [
            ("synthetic-bootstrap-final", Resampling::Bootstrap),
            ("synthetic-small-sample-final", Resampling::SmallSample),
        ] {
            let d = root().join("tests/golden").join(name);
            let r: ReplayManifest =
                serde_json::from_str(&std::fs::read_to_string(d.join("replay.json")).unwrap())
                    .unwrap();
            c.resampling = mode;
            c.bootstrap_iterations = 100;
            let (plans, _) = build_plans(&r.genes, &ess, &non, &c).unwrap();
            assert_eq!(
                serde_json::to_vec(&plans).unwrap(),
                serde_json::to_vec(&r.iterations).unwrap(),
                "{name}"
            );
        }
    }
    #[test]
    #[ignore = "bounded scientific qualification gate, not a routine software test"]
    fn native_bootstrap_distributional_gate() {
        let dir = root().join("tests/external/qualification/bootstrap1000");
        let m = load_matrix(&dir.join("reference.foldchange"));
        let ess = reference(&dir.join("inputs/essentials.tsv"));
        let non = reference(&dir.join("inputs/nonessentials.tsv"));
        let expected = npy(&dir.join("arrays/gene_mean_std_count.npy"));
        let mut c = BfConfig {
            columns: vec![0, 1],
            resampling: Resampling::Bootstrap,
            bootstrap_iterations: 1000,
            rng_mode: RngMode::Rust,
            threads: 4,
            ..Default::default()
        };
        let result = run_bf_base(&m, &ess, &non, &c, None).unwrap();
        let n = result.genes.len();
        let a: Vec<f64> = result.genes.iter().map(|g| g.raw_mean).collect();
        let b: Vec<f64> = expected.as_chunks::<3>().0.iter().map(|v| v[0]).collect();
        let ranks = |v: &[f64]| -> (Vec<f64>, Vec<usize>) {
            let mut order: Vec<usize> = (0..v.len()).collect();
            order.sort_by(|&i, &j| v[j].total_cmp(&v[i]));
            let mut ranks = vec![0.; v.len()];
            let mut i = 0;
            while i < order.len() {
                let mut j = i + 1;
                while j < order.len() && v[order[j]] == v[order[i]] {
                    j += 1;
                }
                let rank = (i + j - 1) as f64 / 2.;
                for k in i..j {
                    ranks[order[k]] = rank;
                }
                i = j;
            }
            (ranks, order)
        };
        let (ra, oa) = ranks(&a);
        let (rb, ob) = ranks(&b);
        let ma = mean(&ra);
        let mb = mean(&rb);
        let cov: f64 = ra.iter().zip(&rb).map(|(x, y)| (x - ma) * (y - mb)).sum();
        let va: f64 = ra.iter().map(|x| (x - ma).powi(2)).sum();
        let vb: f64 = rb.iter().map(|x| (x - mb).powi(2)).sum();
        let spearman = cov / (va * vb).sqrt();
        let overlap = oa[..100].iter().filter(|i| ob[..100].contains(i)).count();
        eprintln!("native bootstrap1000 n={n}: Spearman={spearman:.15}, top100_overlap={overlap}");
        assert!(spearman >= 0.999, "native Spearman {spearman}");
        assert!(overlap >= 98, "native overlap {overlap}");
        // Native CV is deterministic independent of worker count as well.
        c.resampling = Resampling::CrossValidation;
        c.threads = 1;
        let expected = serde_json::to_vec(&run_bf_base(&m, &ess, &non, &c, None).unwrap()).unwrap();
        for t in [2, 4, 8] {
            c.threads = t;
            assert_eq!(
                expected,
                serde_json::to_vec(&run_bf_base(&m, &ess, &non, &c, None).unwrap()).unwrap()
            );
        }
    }
    #[test]
    #[ignore = "explicit HAP1/bootstrap1000 native G7 evidence generation"]
    fn capture_native_hap1_1000() {
        capture_native_1000(false, RngMode::Rust);
    }
    #[test]
    #[ignore = "explicit synthetic/bootstrap1000 native G7 evidence generation"]
    fn capture_native_synthetic_1000() {
        capture_native_1000(true, RngMode::Rust);
    }
    #[test]
    #[ignore = "explicit NumPy-compatible HAP1/bootstrap1000 kernel diagnosis"]
    fn capture_numpy_hap1_1000() {
        capture_native_1000(false, RngMode::Numpy);
    }
    fn capture_native_1000(synthetic: bool, rng_mode: RngMode) {
        use std::io::Write;
        let output = root().join("tests/external/qualification/native-rng");
        let dir = if synthetic {
            root().join("tests/external/qualification/bootstrap1000")
        } else {
            root().join("tests/golden/hap1-cv-final")
        };
        let m = load_matrix(&dir.join("reference.foldchange"));
        let ess = reference(&if synthetic {
            dir.join("inputs/essentials.tsv")
        } else {
            root().join("upstream/BAGEL2/CEGv2.txt")
        });
        let non = reference(&if synthetic {
            dir.join("inputs/nonessentials.tsv")
        } else {
            root().join("upstream/BAGEL2/NEGv1.txt")
        });
        let prefix = if rng_mode == RngMode::Numpy {
            "numpy-compatible"
        } else if synthetic {
            "native-synthetic"
        } else {
            "native"
        };
        let config = BfConfig {
            columns: (0..m.cols()).collect(),
            resampling: Resampling::Bootstrap,
            bootstrap_iterations: 1000,
            rng_mode,
            threads: 4,
            ..Default::default()
        };
        let result = run_bf_base(&m, &ess, &non, &config, None).unwrap();
        let values: Vec<f64> = result
            .genes
            .iter()
            .flat_map(|g| [g.raw_mean, g.std, g.count as f64])
            .collect();
        assert!(values.iter().all(|v| v.is_finite()));
        let header = format!(
            "{{'descr': '<f8', 'fortran_order': False, 'shape': ({}, 3), }}",
            result.genes.len()
        );
        let padding = (64 - ((10 + header.len() + 1) % 64)) % 64;
        let header = header + &" ".repeat(padding) + "\n";
        let mut file = std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(output.join(format!("{prefix}-gene-mean-std-count.npy")))
            .unwrap();
        file.write_all(b"\x93NUMPY\x01\x00").unwrap();
        file.write_all(&(header.len() as u16).to_le_bytes())
            .unwrap();
        file.write_all(header.as_bytes()).unwrap();
        for value in values {
            file.write_all(&value.to_le_bytes()).unwrap();
        }
        let genes: Vec<&str> = result.genes.iter().map(|g| g.gene.as_str()).collect();
        let mut file = std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(output.join(format!("{prefix}-genes.json")))
            .unwrap();
        serde_json::to_writer(&mut file, &genes).unwrap();
        let mut file = std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(output.join(format!("{prefix}-metadata.json")))
            .unwrap();
        serde_json::to_writer_pretty(&mut file,&serde_json::json!({"kind":"scientific-qualification-not-benchmark","config":config,"genes":result.genes.len(),"guides":result.guides.len(),"all_raw_summary_finite":true})).unwrap();
        eprintln!(
            "{prefix}/bootstrap1000 capture complete: {} genes",
            result.genes.len()
        );
    }
    #[test]
    fn cancellation_before_dispatch() {
        let dir = root().join("tests/golden/synthetic-cv-final-a");
        let m = load_matrix(&dir.join("reference.foldchange"));
        let c = BfConfig {
            columns: vec![0, 1],
            ..Default::default()
        };
        assert!(matches!(
            run_bf_base(
                &m,
                &reference(&dir.join("inputs/essentials.tsv")),
                &reference(&dir.join("inputs/nonessentials.tsv")),
                &c,
                Some(&|_, _| ControlFlow::Break(()))
            ),
            Err(Error::Cancelled)
        ));
    }
}
