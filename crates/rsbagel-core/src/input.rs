use crate::bf_types::AlignmentRecord;
use crate::{Error, Matrix, Result};
use std::io::Read;

fn missing(s: &str) -> bool {
    matches!(
        s,
        "" | "#N/A"
            | "#N/A N/A"
            | "#NA"
            | "-1.#IND"
            | "-1.#QNAN"
            | "-NaN"
            | "-nan"
            | "1.#IND"
            | "1.#QNAN"
            | "<NA>"
            | "N/A"
            | "NA"
            | "NULL"
            | "NaN"
            | "n/a"
            | "nan"
            | "null"
    )
}
/// Parse numeric matrices without imposing duplicate-guide policy (BF does so after exclusion).
pub fn read_matrix(reader: impl Read, counts: bool) -> Result<Matrix> {
    let _profile_span = crate::profiling::span("parse_matrix");
    let mut csv = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .quoting(counts)
        .flexible(false)
        .from_reader(reader);
    let header = csv.headers()?.clone();
    if header.len() < 3 {
        return Err(Error::Invalid(
            "matrix requires reagent, gene and numeric columns".into(),
        ));
    }
    // pandas mangles duplicate headers; BF's manual reader retains them.
    let mut columns = Vec::new();
    for label in header.iter().skip(2) {
        let mut name = label.to_string();
        if counts {
            let mut suffix = 1;
            while columns.contains(&name) {
                name = format!("{label}.{suffix}");
                suffix += 1;
            }
        }
        columns.push(name);
    }
    let mut m = Matrix {
        index_name: header[0].to_string(),
        row_ids: vec![],
        row_genes: vec![],
        columns,
        values: vec![],
    };
    for record in csv.records() {
        let record = record?;
        m.row_ids.push(record[0].to_string());
        m.row_genes.push(if counts && missing(&record[1]) {
            "NO_GENE_NAME".into()
        } else {
            record[1].to_string()
        });
        for field in record.iter().skip(2) {
            if counts && missing(field) {
                m.values.push(0.0);
            } else {
                m.values.push(field.trim().parse::<f64>().map_err(|_| {
                    Error::Invalid(format!("could not convert string to float: '{field}'"))
                })?);
            }
        }
    }
    m.validate()?;
    Ok(m)
}
pub fn read_references(mut reader: impl Read) -> Result<Vec<String>> {
    let _profile_span = crate::profiling::span("parse_references");
    let mut text = String::new();
    reader.read_to_string(&mut text)?;
    Ok(text
        .lines()
        .skip(1)
        .map(|s| s.trim_end().split('\t').next().unwrap_or("").to_string())
        .collect())
}
pub fn read_alignment(mut reader: impl Read) -> Result<Vec<AlignmentRecord>> {
    let mut text = String::new();
    reader.read_to_string(&mut text)?;
    let mut result = vec![];
    for line in text.lines() {
        let fields: Vec<_> = line.split('\t').collect();
        if fields.len() != 5 {
            return Err(Error::Invalid("Please check align-info file".into()));
        }
        let count = |s: &str| if missing(s) { 0 } else { s.split(',').count() };
        result.push(AlignmentRecord {
            guide: fields[0].into(),
            perfect_loci: count(fields[1]),
            perfect_genes: count(fields[2]),
            mismatch_loci: count(fields[3]),
            mismatch_genes: count(fields[4]),
        });
    }
    let mut ids = std::collections::HashSet::new();
    if result.iter().any(|x| !ids.insert(x.guide.clone())) {
        return Err(Error::Invalid("Please check align-info file".into()));
    }
    Ok(result)
}
pub fn read_network(mut reader: impl Read) -> Result<Vec<(String, String)>> {
    let _profile_span = crate::profiling::span("parse_network");
    let mut text = String::new();
    reader.read_to_string(&mut text)?;
    text.lines()
        .map(|line| {
            let f: Vec<_> = line.trim_end().split('\t').collect();
            if f.len() < 2 {
                Err(Error::Invalid("network row requires two endpoints".into()))
            } else {
                Ok((f[0].into(), f[1].into()))
            }
        })
        .collect()
}
/// Upstream numbers index columns after the guide index: gene=0, first numeric=1.
pub fn select_columns(m: &Matrix, selector: &str, bf: bool) -> Result<(Vec<usize>, Vec<String>)> {
    let tokens: Vec<_> = selector.split(',').collect();
    if let Ok(numbers) = tokens
        .iter()
        .map(|x| x.parse::<isize>())
        .collect::<std::result::Result<Vec<_>, _>>()
    {
        let mut out = vec![];
        let mut labels = vec![];
        for x in numbers {
            let actual = if x < 0 {
                m.cols() as isize + if bf { 2 } else { 1 } + x
            } else {
                x
            };
            if (bf && x == -1) || actual < 1 || actual > m.cols() as isize {
                return Err(Error::Invalid("Invalid columns".into()));
            }
            let i = (actual - 1) as usize;
            out.push(i);
            labels.push(m.columns[i].clone());
        }
        Ok((out, labels))
    } else {
        let labels = tokens.iter().map(|s| s.to_string()).collect::<Vec<_>>();
        if bf {
            Ok((
                m.columns
                    .iter()
                    .enumerate()
                    .filter(|(_, x)| labels.contains(x))
                    .map(|(i, _)| i)
                    .collect(),
                labels,
            ))
        } else {
            let columns = labels
                .iter()
                .map(|s| {
                    m.columns
                        .iter()
                        .position(|x| x == s)
                        .ok_or_else(|| Error::Invalid("Invalid input controls".into()))
                })
                .collect::<Result<Vec<_>>>()?;
            Ok((columns, labels))
        }
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn pinned_parser_label_and_selector_regressions() {
        let matrix = super::read_matrix(
            &include_bytes!("../../../tests/fixtures/parser/counts.tsv")[..],
            true,
        )
        .unwrap();
        assert_eq!(matrix.row_genes, ["None", "NO_GENE_NAME", "A"]);
        assert_eq!(super::select_columns(&matrix, "-1", false).unwrap().0, [1]);
        assert_eq!(super::select_columns(&matrix, "-2", true).unwrap().0, [1]);
        assert!(super::select_columns(&matrix, "-1", true).is_err());
        let quoted = &include_bytes!("../../../tests/fixtures/parser/quoted.tsv")[..];
        assert!(super::read_matrix(quoted, false).is_err());
    }
}
