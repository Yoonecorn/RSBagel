use crate::{Error, Matrix, Result};
use serde::{Deserialize, Serialize};
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FcConfig {
    pub controls: Vec<usize>,
    pub control_labels: Vec<String>,
    pub min_reads: f64,
    pub pseudocount: f64,
}
impl Default for FcConfig {
    fn default() -> Self {
        Self {
            controls: vec![0],
            control_labels: vec![],
            min_reads: 0.0,
            pseudocount: 5.0,
        }
    }
}
#[derive(Clone, Debug, Serialize)]
pub struct FcResult {
    pub fold_change: Matrix,
    pub normalized: Matrix,
}
pub fn run_fc(m: &Matrix, cfg: &FcConfig) -> Result<FcResult> {
    let _profile_span = crate::profiling::span("fold_change");
    m.validate()?;
    if [cfg.min_reads, cfg.pseudocount]
        .into_iter()
        .any(|v| !v.is_finite() || v.fract() != 0.0)
    {
        return Err(Error::Invalid(
            "minimum reads and pseudocount must be finite integers".into(),
        ));
    }
    if cfg.controls.is_empty() || cfg.controls.iter().any(|&x| x >= m.cols()) {
        return Err(Error::Invalid("Invalid input controls".into()));
    }
    let labels = if cfg.control_labels.is_empty() {
        cfg.controls
            .iter()
            .map(|&x| m.columns[x].clone())
            .collect::<Vec<_>>()
    } else {
        cfg.control_labels.clone()
    };
    let control_label = labels.join(";");
    let remaining = (0..m.cols())
        .filter(|i| !cfg.controls.contains(i))
        .collect::<Vec<_>>();
    let mut columns = remaining
        .iter()
        .map(|&i| m.columns[i].clone())
        .collect::<Vec<_>>();
    let control_pos = if let Some(i) = columns.iter().position(|s| s == &control_label) {
        i
    } else {
        columns.push(control_label.clone());
        columns.len() - 1
    };
    let c = columns.len();
    let mut added = vec![0.0; m.rows() * c];
    let mut totals = vec![0.0; c];
    for r in 0..m.rows() {
        for (j, &i) in remaining.iter().enumerate() {
            added[r * c + j] = m.at(r, i) + cfg.pseudocount;
        }
        added[r * c + control_pos] =
            cfg.controls.iter().map(|&i| m.at(r, i)).sum::<f64>() + cfg.pseudocount;
        for j in 0..c {
            totals[j] += added[r * c + j];
        }
    }
    let mut normalized = Matrix {
        index_name: m.index_name.clone(),
        row_ids: vec![],
        row_genes: vec![],
        columns: columns.clone(),
        values: vec![],
    };
    let mut fc = Matrix {
        index_name: "REAGENT_ID".into(),
        row_ids: vec![],
        row_genes: vec![],
        columns: columns
            .iter()
            .enumerate()
            .filter(|(i, _)| *i != control_pos)
            .map(|(_, x)| x.clone())
            .collect(),
        values: vec![],
    };
    for r in 0..m.rows() {
        if added[r * c + control_pos] < cfg.min_reads || added[r * c + control_pos].is_nan() {
            continue;
        }
        normalized.row_ids.push(m.row_ids[r].clone());
        normalized.row_genes.push(m.row_genes[r].clone());
        fc.row_ids.push(m.row_ids[r].clone());
        fc.row_genes.push(m.row_genes[r].clone());
        let norm = (0..c)
            .map(|j| added[r * c + j] / totals[j] * 10_000_000.0)
            .collect::<Vec<_>>();
        for j in 0..c {
            if j != control_pos {
                fc.values.push((norm[j] / norm[control_pos]).log2());
            }
        }
        normalized.values.extend(norm);
    }
    Ok(FcResult {
        fold_change: fc,
        normalized,
    })
}

/// Match the FC writer -> BF reader precision boundary without disk I/O.
/// Raw FC remains available independently for kernel-level qualification.
pub fn fc_to_bf_input(fold_change: &Matrix) -> Result<Matrix> {
    let _profile_span = crate::profiling::span("memory_fc_precision_bridge");
    fold_change.validate()?;
    let mut input = fold_change.clone();
    for value in &mut input.values {
        if value.is_nan() {
            return Err(Error::Invalid(
                "FC serialization contains an empty numeric field; BF input is invalid".into(),
            ));
        }
        *value = crate::output::float(*value, 3, 4)
            .trim()
            .parse()
            .map_err(|_| Error::Invalid("FC serialization is not a valid BF number".into()))?;
    }
    Ok(input)
}
#[cfg(test)]
mod tests {
    #[test]
    fn memory_precision_bridge_matches_source_bf_input() {
        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/golden/synthetic-cv-final-a");
        let counts = crate::input::read_matrix(
            std::fs::File::open(root.join("inputs/counts.tsv")).unwrap(),
            true,
        )
        .unwrap();
        let result = super::run_fc(&counts, &super::FcConfig::default()).unwrap();
        let bridge = super::fc_to_bf_input(&result.fold_change).unwrap();
        let expected = crate::input::read_matrix(
            std::fs::File::open(root.join("reference.foldchange")).unwrap(),
            false,
        )
        .unwrap();
        assert_eq!(bridge.row_ids, expected.row_ids);
        assert_eq!(bridge.row_genes, expected.row_genes);
        assert_eq!(bridge.columns, expected.columns);
        assert_eq!(
            bridge
                .values
                .iter()
                .map(|x| x.to_bits())
                .collect::<Vec<_>>(),
            expected
                .values
                .iter()
                .map(|x| x.to_bits())
                .collect::<Vec<_>>()
        );
        let mut invalid = bridge;
        invalid.values[0] = f64::NAN;
        assert!(super::fc_to_bf_input(&invalid).is_err());
    }
}
