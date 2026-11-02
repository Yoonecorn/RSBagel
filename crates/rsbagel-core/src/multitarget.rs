//! BAGEL2 alignment correction and guide-count equalization.
//! Copyright (c) 2020 Hart Lab; MIT, see docs/BAGEL2-LICENSE.txt.
use crate::bf_types::{AlignmentRecord, BfConfig, BfResult, MultiTargetTrace};
use crate::kde::mean;
use crate::{Error, Result};
use std::collections::HashMap;

fn is_multiple(a: &AlignmentRecord) -> bool {
    a.perfect_loci > 1 || a.mismatch_loci > 0
}

/// Least squares with an intercept, using a two-column Jacobi SVD. The intercept
/// is fitted through centering but, as upstream, is not applied in the penalty.
fn fit_two_predictors(x: &[[f64; 2]], y: &[f64]) -> Result<[f64; 2]> {
    if x.is_empty() || x.len() != y.len() || y.iter().any(|v| !v.is_finite()) {
        return Err(Error::Numerical(
            "invalid multi-target regression observations".into(),
        ));
    }
    let first: Vec<_> = x.iter().map(|v| v[0]).collect();
    let second: Vec<_> = x.iter().map(|v| v[1]).collect();
    let mx = [mean(&first), mean(&second)];
    let my = mean(y);
    let centered: Vec<_> = x.iter().map(|v| [v[0] - mx[0], v[1] - mx[1]]).collect();
    let mut a = 0.0;
    let mut b = 0.0;
    let mut c = 0.0;
    for v in &centered {
        a += v[0] * v[0];
        b += v[0] * v[1];
        c += v[1] * v[1];
    }
    let angle = 0.5 * (2.0 * b).atan2(a - c);
    let (sn, cs) = angle.sin_cos();
    let vectors = [[cs, sn], [-sn, cs]];
    let mut singular2 = [0.0; 2];
    let mut rhs = [0.0; 2];
    for (row, &response) in centered.iter().zip(y) {
        for j in 0..2 {
            let projected = row[0] * vectors[j][0] + row[1] * vectors[j][1];
            singular2[j] += projected * projected;
            rhs[j] += projected * (response - my);
        }
    }
    let cutoff = f64::EPSILON * singular2[0].max(singular2[1]).sqrt();
    let mut coefficients = [0.0; 2];
    for j in 0..2 {
        if singular2[j].sqrt() > cutoff {
            let coefficient = rhs[j] / singular2[j];
            coefficients[0] += coefficient * vectors[j][0];
            coefficients[1] += coefficient * vectors[j][1];
        }
    }
    if coefficients.iter().any(|v| !v.is_finite()) {
        return Err(Error::Numerical(
            "non-finite multi-target regression coefficients".into(),
        ));
    }
    Ok(coefficients)
}

/// Applies only the upstream score correction/equalization. Physical alignment
/// exclusion is already performed before duplicate guide validation in BF base.
pub fn apply_multitarget(result: &mut BfResult, config: &BfConfig) -> Result<()> {
    if let Some(alignment) = &config.alignment {
        let mut by_guide = HashMap::new();
        for a in alignment {
            if by_guide.insert(a.guide.as_str(), a).is_some() {
                return Err(Error::Legacy("Please check align-info file".into()));
            }
        }
        let mut x = Vec::new();
        let mut y = Vec::new();
        let mut eligible = 0;
        let mut training_guides = Vec::new();
        for indices in &result.gene_guides {
            let mut only = Vec::new();
            let mut multiple = Vec::new();
            for &index in indices {
                if let Some(a) = by_guide.get(result.guides[index].guide.as_str()) {
                    if is_multiple(a) {
                        multiple.push((index, *a));
                    } else {
                        only.push(index);
                    }
                }
            }
            if only.is_empty() || multiple.is_empty() {
                continue;
            }
            let singles: Vec<_> = only.iter().map(|&i| result.guides[i].total).collect();
            let baseline = mean(&singles);
            for index in only {
                training_guides.push(result.guides[index].guide.clone());
                x.push([1.0, 0.0]);
                y.push(0.0);
            }
            for (index, a) in multiple {
                if a.perfect_genes > 1 || a.mismatch_genes > 0 {
                    continue;
                }
                eligible += 1;
                training_guides.push(result.guides[index].guide.clone());
                x.push([a.perfect_loci as f64, a.mismatch_loci as f64]);
                y.push(result.guides[index].total - baseline);
            }
        }
        if eligible < 10 {
            result.diagnostics.extend([
                "Not enough train set for calculating multi-targeting effect.\n".into(),
                "It may cause due to unmatched gRNA names between the foldchange file and the align info file.\n".into(),
                "Filtering is not finished\n".into(),
            ]);
        } else {
            let mut coefficients = fit_two_predictors(&x, &y)?;
            let raw_coefficients = coefficients;
            for coefficient in &mut coefficients {
                if *coefficient < 0.0 {
                    result.diagnostics.push(
                        "Regression coefficient is below than zero. Substituted to zero\n".into(),
                    );
                    *coefficient = 0.0;
                }
            }
            result.diagnostics.push(format!("Multiple effects from perfect matched loci = {:.3} and 1bp mis-matched loci = {:.3}", coefficients[0], coefficients[1]));
            result.correction_coefficients = Some(coefficients);
            if config.dump_stages {
                let mx0 = mean(&x.iter().map(|v| v[0]).collect::<Vec<_>>());
                let mx1 = mean(&x.iter().map(|v| v[1]).collect::<Vec<_>>());
                let intercept = mean(&y) - raw_coefficients[0] * mx0 - raw_coefficients[1] * mx1;
                result.correction_trace = Some(MultiTargetTrace {
                    guides: training_guides,
                    x,
                    y,
                    coefficients_unclipped: raw_coefficients,
                    coefficients_clipped: coefficients,
                    intercept,
                });
            }
            for (g, indices) in result.gene_guides.iter().enumerate() {
                let mut gene_penalty = 0.0;
                for &index in indices {
                    let penalty = by_guide
                        .get(result.guides[index].guide.as_str())
                        .filter(|a| is_multiple(a))
                        .map_or(0.0, |a| {
                            (a.perfect_loci as f64 - 1.0) * coefficients[0]
                                + a.mismatch_loci as f64 * coefficients[1]
                        });
                    gene_penalty += penalty;
                    if config.guide_level {
                        result.guides[index].total -= penalty;
                    }
                }
                if !config.guide_level {
                    result.genes[g].mean = result.genes[g].raw_mean - gene_penalty;
                }
            }
        }
    }
    if !config.guide_level {
        if let Some(count) = config.equalise_guides.filter(|v| *v != 0.0) {
            for (gene, indices) in result.genes.iter_mut().zip(&result.gene_guides) {
                gene.normalized = Some(gene.mean * (count / indices.len() as f64));
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn two_predictor_ols_and_rank_deficiency() {
        let x = [[1., 0.], [2., 0.], [1., 1.], [3., 2.], [4., 1.]];
        let y: Vec<_> = x.iter().map(|v| 7. + 2. * v[0] - 3. * v[1]).collect();
        let b = fit_two_predictors(&x, &y).unwrap();
        assert!((b[0] - 2.).abs() < 1e-12 && (b[1] + 3.).abs() < 1e-12);
        let x = [[1., 0.], [2., 0.], [3., 0.]];
        let b = fit_two_predictors(&x, &[3., 5., 7.]).unwrap();
        assert!((b[0] - 2.).abs() < 1e-12 && b[1].abs() < 1e-12);
        let x = [[1., 2.], [2., 4.], [3., 6.]];
        let b = fit_two_predictors(&x, &[3., 6., 9.]).unwrap();
        assert!((b[0] - 0.6).abs() < 1e-12 && (b[1] - 1.2).abs() < 1e-12);
    }
    #[test]
    fn positive_penalty_preserves_raw_values_and_separate_normalization() {
        use crate::bf_types::{GeneScore, GuideScore};
        let mut genes = Vec::new();
        let mut guides = Vec::new();
        let mut gene_guides = Vec::new();
        let mut alignment = Vec::new();
        for i in 0..12 {
            let name = format!("G{i:03}");
            let base = i as f64 / 3.;
            let raw = 3. * base + 1.6;
            genes.push(GeneScore {
                gene: name.clone(),
                raw_mean: raw,
                mean: raw,
                std: 0.7,
                count: 1,
                normalized: None,
                boosted: None,
                boosted_std: None,
                boosted_normalized: None,
                network_evidence: None,
                network_contribution_count: None,
                network_zero_applied: None,
            });
            let mut indices = Vec::new();
            for (j, (perfect, mismatch, total)) in
                [(1, 0, base), (2, 0, base + 0.4), (1, 1, base + 1.2)]
                    .into_iter()
                    .enumerate()
            {
                let guide = format!("{name}_g{j}");
                indices.push(guides.len());
                guides.push(GuideScore {
                    guide: guide.clone(),
                    gene: name.clone(),
                    mean_replicates: vec![total],
                    std_replicates: vec![0.],
                    total,
                    count: 1,
                });
                alignment.push(AlignmentRecord {
                    guide,
                    perfect_loci: perfect,
                    perfect_genes: 1,
                    mismatch_loci: mismatch,
                    mismatch_genes: 0,
                });
            }
            gene_guides.push(indices);
        }
        let mut result = BfResult {
            genes,
            guides,
            gene_guides,
            partitions: Vec::new(),
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
        let config = BfConfig {
            alignment: Some(alignment),
            equalise_guides: Some(4.),
            ..BfConfig::default()
        };
        apply_multitarget(&mut result, &config).unwrap();
        let coefficients = result.correction_coefficients.unwrap();
        assert!((coefficients[0] - 0.4).abs() < 1e-12);
        assert!((coefficients[1] - 1.2).abs() < 1e-12);
        for (i, gene) in result.genes.iter().enumerate() {
            assert!((gene.mean - i as f64).abs() < 1e-12);
            assert!((gene.raw_mean - gene.mean - 1.6).abs() < 1e-12);
            assert_eq!(gene.std, 0.7);
            assert!((gene.normalized.unwrap() - gene.mean * 4. / 3.).abs() < 1e-12);
        }
    }
}
