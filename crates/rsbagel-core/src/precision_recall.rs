use crate::{Error, Result};
use serde::Serialize;
use std::collections::HashSet;
#[derive(Clone, Debug, Serialize)]
pub struct PrRow {
    pub gene: String,
    pub score: f64,
    pub recall: f64,
    pub precision: f64,
    pub fdr: f64,
}

// Port of NumPy1.23.5 numeric indirect introsort. Ties follow the pinned
// unstable permutation rather than an invented stable/alphabetical policy.
fn heap(a: &mut [usize], v: &[f64]) {
    let n = a.len();
    if n < 2 {
        return;
    }
    let mut l = n / 2;
    let mut ir = n;
    loop {
        let tmp = if l > 0 {
            l -= 1;
            a[l]
        } else {
            ir -= 1;
            let t = a[ir];
            a[ir] = a[0];
            if ir == 1 {
                a[0] = t;
                break;
            }
            t
        };
        let mut i = l;
        let mut j = 2 * l + 1;
        while j < ir {
            if j + 1 < ir && v[a[j]] < v[a[j + 1]] {
                j += 1;
            }
            if v[tmp] < v[a[j]] {
                a[i] = a[j];
                i = j;
                j = 2 * j + 1;
            } else {
                break;
            }
        }
        a[i] = tmp;
    }
}
fn argsort(v: &[f64]) -> Vec<usize> {
    let mut a = (0..v.len()).collect::<Vec<_>>();
    if a.len() < 2 {
        return a;
    }
    let mut stack = vec![(
        0isize,
        a.len() as isize - 1,
        2 * (usize::BITS - a.len().leading_zeros() - 1) as i32,
    )];
    while let Some((mut lo, mut hi, mut depth)) = stack.pop() {
        if depth < 0 {
            if hi >= lo {
                heap(&mut a[lo as usize..=hi as usize], v);
            }
            continue;
        }
        while hi - lo > 15 {
            let mid = lo + ((hi - lo) >> 1);
            if v[a[mid as usize]] < v[a[lo as usize]] {
                a.swap(mid as usize, lo as usize);
            }
            if v[a[hi as usize]] < v[a[mid as usize]] {
                a.swap(hi as usize, mid as usize);
            }
            if v[a[mid as usize]] < v[a[lo as usize]] {
                a.swap(mid as usize, lo as usize);
            }
            let pivot = v[a[mid as usize]];
            let mut i = lo;
            let mut j = hi - 1;
            a.swap(mid as usize, j as usize);
            loop {
                i += 1;
                while v[a[i as usize]] < pivot {
                    i += 1;
                }
                j -= 1;
                while pivot < v[a[j as usize]] {
                    j -= 1;
                }
                if i >= j {
                    break;
                }
                a.swap(i as usize, j as usize);
            }
            a.swap(i as usize, (hi - 1) as usize);
            depth -= 1;
            if i - lo < hi - i {
                stack.push((i + 1, hi, depth));
                hi = i - 1;
            } else {
                stack.push((lo, i - 1, depth));
                lo = i + 1;
            }
        }
        for i in (lo + 1)..=hi {
            let item = a[i as usize];
            let mut j = i;
            while j > lo && v[item] < v[a[(j - 1) as usize]] {
                a[j as usize] = a[(j - 1) as usize];
                j -= 1;
            }
            a[j as usize] = item;
        }
    }
    a
}
/// pandas1.5.3 descending nargsort: reverse non-NaN inputs, argsort, reverse indices.
pub fn descending_order(scores: &[f64]) -> Vec<usize> {
    let mut indices = (0..scores.len())
        .filter(|&i| !scores[i].is_nan())
        .collect::<Vec<_>>();
    indices.reverse();
    let vals = indices.iter().map(|&i| scores[i]).collect::<Vec<_>>();
    let mut order = argsort(&vals)
        .iter()
        .map(|&i| indices[i])
        .collect::<Vec<_>>();
    order.reverse();
    order.extend((0..scores.len()).filter(|&i| scores[i].is_nan()));
    order
}
pub fn run_pr(
    genes: &[String],
    scores: &[f64],
    ess: &[String],
    non: &[String],
) -> Result<Vec<PrRow>> {
    let _profile_span = crate::profiling::span("precision_recall");
    if genes.len() != scores.len() {
        return Err(Error::Invalid("gene and score lengths differ".into()));
    }
    let mut ids = HashSet::new();
    if genes.iter().any(|g| !ids.insert(g)) {
        return Err(Error::Legacy(
            "unsupported format string passed to Series.__format__".into(),
        ));
    }
    let ess: HashSet<_> = ess.iter().collect();
    let non: HashSet<_> = non.iter().collect();
    let denominator = genes.iter().filter(|g| ess.contains(g)).count() as f64;
    if denominator == 0.0 && !genes.is_empty() {
        return Err(Error::Legacy("float division by zero".into()));
    }
    let mut tp = 0.0;
    let mut fp = 0.0;
    let mut precision = 1.0;
    let mut result = vec![];
    for i in descending_order(scores) {
        let g = &genes[i];
        if ess.contains(g) {
            tp += 1.0;
        } else if non.contains(g) {
            fp += 1.0;
        }
        if tp > 0.0 || fp > 0.0 {
            precision = tp / (tp + fp);
        }
        result.push(PrRow {
            gene: g.clone(),
            score: scores[i],
            recall: tp / denominator,
            precision,
            fdr: 1.0 - precision,
        });
    }
    Ok(result)
}

/// Reproduce the alphabetical BF writer and its three-decimal numeric boundary
/// before PR, without disk serialization. Column selection matches CLI `pr -k`.
pub fn bf_to_pr_input(
    result: &crate::BfResult,
    cfg: &crate::BfConfig,
    column: &str,
) -> Result<(Vec<String>, Vec<f64>)> {
    let _profile_span = crate::profiling::span("memory_bf_precision_bridge");
    use crate::{NetworkProfile, Resampling};
    if cfg.guide_level {
        return Err(Error::Invalid(
            "gene PR requires gene-level BF output".into(),
        ));
    }
    let bootstrap = cfg.resampling != Resampling::CrossValidation;
    let network = cfg.network.is_some();
    let eqf = crate::output::replicate_factor(cfg);
    let mut rows = result.genes.iter().collect::<Vec<_>>();
    rows.sort_by(|a, b| a.gene.cmp(&b.gene));
    let mut genes = Vec::with_capacity(rows.len());
    let mut scores = Vec::with_capacity(rows.len());
    for g in rows {
        let value = match column {
            "BF" => g.mean * eqf,
            "STD" if bootstrap => g.std,
            "NumObs" if bootstrap => g.count as f64,
            "NormBF" if cfg.equalise_guides.is_some_and(|v| v != 0.0) => g
                .normalized
                .ok_or_else(|| Error::Invalid("NormBF unavailable".into()))?,
            "BoostedBF" if network => {
                g.boosted
                    .ok_or_else(|| Error::Invalid("BoostedBF unavailable".into()))?
                    * eqf
            }
            "STD_BoostedBF" if network && bootstrap => {
                g.boosted_std
                    .ok_or_else(|| Error::Invalid("STD_BoostedBF unavailable".into()))?
                    * eqf
            }
            "BoostedNormBF"
                if network
                    && cfg.network_profile == NetworkProfile::Corrected
                    && cfg.equalise_guides.is_some_and(|v| v != 0.0) =>
            {
                g.boosted_normalized
                    .ok_or_else(|| Error::Invalid("BoostedNormBF unavailable".into()))?
                    * eqf
            }
            _ => {
                return Err(Error::Invalid(
                    "Error! the column name is not in the file".into(),
                ))
            }
        };
        let score = if value.is_nan() {
            f64::NAN
        } else {
            crate::output::float(value, 3, 4)
                .trim()
                .parse()
                .map_err(|_| Error::Invalid("BF serialization is not a valid PR number".into()))?
        };
        genes.push(g.gene.clone());
        scores.push(score);
    }
    Ok((genes, scores))
}
#[cfg(test)]
mod tests {
    #[test]
    fn memory_bf_pr_bridge_matches_source_output_and_order() {
        use crate::{BfConfig, RngMode};
        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/golden/synthetic-cv-final-a");
        let matrix = crate::input::read_matrix(
            std::fs::File::open(root.join("reference.foldchange")).unwrap(),
            false,
        )
        .unwrap();
        let essentials = crate::input::read_references(
            std::fs::File::open(root.join("inputs/essentials.tsv")).unwrap(),
        )
        .unwrap();
        let non = crate::input::read_references(
            std::fs::File::open(root.join("inputs/nonessentials.tsv")).unwrap(),
        )
        .unwrap();
        let cfg = BfConfig {
            columns: vec![0, 1],
            column_labels: matrix.columns.clone(),
            rng_mode: RngMode::Numpy,
            ..BfConfig::default()
        };
        let mut result = crate::run_bf(&matrix, &essentials, &non, &cfg, None).unwrap();
        result.genes.reverse();
        let (genes, scores) = super::bf_to_pr_input(&result, &cfg, "BF").unwrap();
        let source = std::fs::read_to_string(root.join("reference.bf")).unwrap();
        let expected = source
            .lines()
            .skip(1)
            .map(|line| {
                let mut f = line.split('\t');
                (
                    f.next().unwrap().to_string(),
                    f.next().unwrap().parse::<f64>().unwrap(),
                )
            })
            .collect::<Vec<_>>();
        assert_eq!(
            genes,
            expected.iter().map(|x| x.0.clone()).collect::<Vec<_>>()
        );
        assert_eq!(
            scores.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
            expected.iter().map(|x| x.1.to_bits()).collect::<Vec<_>>()
        );
        let rows = super::run_pr(&genes, &scores, &essentials, &non).unwrap();
        let mut text = Vec::new();
        crate::output::write_pr(&mut text, &rows, "BF").unwrap();
        assert_eq!(text, std::fs::read(root.join("reference.pr")).unwrap());
    }
}
