use crate::bf_types::*;
use crate::precision_recall::PrRow;
use crate::{Matrix, Result};
use std::io::Write;
/// Python fixed-point formatting; preserve lowercase nonfinite values and field widths.
pub fn float(v: f64, precision: usize, width: usize) -> String {
    let text = if v.is_nan() {
        "nan".into()
    } else if v == f64::INFINITY {
        "inf".into()
    } else if v == f64::NEG_INFINITY {
        "-inf".into()
    } else {
        format!("{v:.precision$}")
    };
    format!("{text:>width$}")
}
pub fn write_matrix(w: impl Write, m: &Matrix, normalized: bool) -> Result<()> {
    let _profile_span = crate::profiling::span("format_write_matrix");
    let mut w = csv::WriterBuilder::new().delimiter(b'\t').from_writer(w);
    let mut header = if normalized {
        vec![m.index_name.clone()]
    } else {
        vec!["REAGENT_ID".into(), "GENE".into()]
    };
    header.extend(m.columns.clone());
    w.write_record(header)?;
    for r in 0..m.rows() {
        let mut row = vec![m.row_ids[r].clone()];
        if !normalized {
            row.push(m.row_genes[r].clone());
        }
        for c in 0..m.cols() {
            let v = m.at(r, c);
            row.push(if v.is_nan() {
                String::new()
            } else {
                float(
                    v,
                    if normalized { 2 } else { 3 },
                    if normalized { 3 } else { 4 },
                )
            });
        }
        w.write_record(row)?;
    }
    w.flush()?;
    Ok(())
}
pub fn default_pr_column(cfg: &BfConfig) -> &'static str {
    if cfg.network.is_some() && !cfg.guide_level && cfg.network_profile == NetworkProfile::Corrected
    {
        if cfg.equalise_guides.is_some_and(|v| v != 0.0) {
            "BoostedNormBF"
        } else {
            "BoostedBF"
        }
    } else {
        "BF"
    }
}
pub fn replicate_factor(cfg: &BfConfig) -> f64 {
    let count = if cfg.column_labels.is_empty() {
        cfg.columns.len()
    } else {
        cfg.column_labels.len()
    };
    cfg.equalise_reps
        .filter(|v| *v != 0.0)
        .map(|v| v / count as f64)
        .unwrap_or(1.0)
}
pub fn write_bf(mut w: impl Write, result: &BfResult, cfg: &BfConfig) -> Result<()> {
    let _profile_span = crate::profiling::span("format_write_bf");
    let boot = cfg.resampling != Resampling::CrossValidation;
    let network = cfg.network.is_some() && !cfg.guide_level;
    let eqf = replicate_factor(cfg);
    let norm = cfg.equalise_guides.is_some_and(|x| x != 0.0);
    if cfg.guide_level {
        write!(w, "RNA\tGENE")?;
        for label in cfg.column_labels.iter().take(cfg.columns.len()) {
            write!(w, "\t{label}")?;
            if boot {
                write!(w, "\t{label}_STD")?;
            }
        }
        write!(w, "\tBF")?;
        if boot {
            write!(w, "\tNumObs")?;
        }
        writeln!(w)?;
        let mut guides = result.guides.iter().collect::<Vec<_>>();
        guides.sort_by(|a, b| a.guide.cmp(&b.guide));
        for guide in guides {
            write!(w, "{}\t{}", guide.guide, guide.gene)?;
            for (mean, std) in guide.mean_replicates.iter().zip(&guide.std_replicates) {
                write!(w, "\t{}", float(*mean, 3, 4))?;
                if boot {
                    write!(w, "\t{}", float(*std, 3, 4))?;
                }
            }
            write!(w, "\t{}", float(guide.total * eqf, 3, 4))?;
            if boot {
                write!(w, "\t{}", guide.count)?;
            }
            writeln!(w)?;
        }
    } else {
        write!(w, "GENE")?;
        if network {
            write!(w, "\tBoostedBF")?;
            if boot {
                write!(w, "\tSTD_BoostedBF")?;
            }
        }
        write!(w, "\tBF")?;
        if boot {
            write!(w, "\tSTD\tNumObs")?;
        }
        if norm {
            write!(w, "\tNormBF")?;
        }
        if network && norm && cfg.network_profile == NetworkProfile::Corrected {
            write!(w, "\tBoostedNormBF")?;
        }
        writeln!(w)?;
        let mut genes = result.genes.iter().collect::<Vec<_>>();
        genes.sort_by(|a, b| a.gene.cmp(&b.gene));
        for g in genes {
            write!(w, "{}", g.gene)?;
            if network {
                write!(w, "\t{}", float(g.boosted.unwrap_or(f64::NAN) * eqf, 3, 4))?;
                if boot {
                    write!(
                        w,
                        "\t{}",
                        float(g.boosted_std.unwrap_or(f64::NAN) * eqf, 3, 4)
                    )?;
                }
            }
            write!(w, "\t{}", float(g.mean * eqf, 3, 4))?;
            if boot {
                write!(w, "\t{}\t{}", float(g.std, 3, 4), g.count)?;
            }
            if norm {
                write!(w, "\t{}", float(g.normalized.unwrap_or(f64::NAN), 3, 4))?;
            }
            if network && norm && cfg.network_profile == NetworkProfile::Corrected {
                write!(
                    w,
                    "\t{}",
                    float(g.boosted_normalized.unwrap_or(f64::NAN) * eqf, 3, 4)
                )?;
            }
            writeln!(w)?;
        }
    }
    Ok(())
}
pub fn write_pr(mut w: impl Write, rows: &[PrRow], column: &str) -> Result<()> {
    let _profile_span = crate::profiling::span("format_write_pr");
    writeln!(w, "Gene\t{column}\tRecall\tPrecision\tFDR")?;
    for r in rows {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}",
            r.gene,
            float(r.score, 3, 4),
            float(r.recall, 3, 4),
            float(r.precision, 3, 4),
            float(r.fdr, 3, 4)
        )?;
    }
    Ok(())
}
pub fn write_traininfo(mut w: impl Write, result: &BfResult) -> Result<()> {
    writeln!(w, "#1: Loopcount\n#2: Training set\n#3: Testset")?;
    for p in &result.partitions {
        let names = |ix: &[usize]| {
            ix.iter()
                .map(|&i| result.genes[i].gene.as_str())
                .collect::<Vec<_>>()
                .join(",")
        };
        writeln!(
            w,
            "{}\n{}\n{}",
            p.iteration,
            names(&p.train),
            names(&p.test)
        )?;
    }
    Ok(())
}
pub fn write_netscore(mut w: impl Write, result: &BfResult) -> Result<()> {
    for e in &result.network_events {
        writeln!(
            w,
            "{}\t{}\t{}",
            e.gene,
            float(e.score, 6, 0),
            float(e.contribution, 6, 0)
        )?;
    }
    Ok(())
}
