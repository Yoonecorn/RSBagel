use serde::{Deserialize, Serialize};

/// Row-major f64 matrix; row identifiers and gene labels retain source order.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Matrix {
    #[serde(default)]
    pub index_name: String,
    pub row_ids: Vec<String>,
    pub row_genes: Vec<String>,
    pub columns: Vec<String>,
    pub values: Vec<f64>,
}
impl Matrix {
    pub fn rows(&self) -> usize {
        self.row_ids.len()
    }
    pub fn cols(&self) -> usize {
        self.columns.len()
    }
    pub fn at(&self, row: usize, column: usize) -> f64 {
        self.values[row * self.cols() + column]
    }
    pub fn validate(&self) -> crate::Result<()> {
        if self.row_genes.len() != self.rows()
            || self.values.len()
                != self
                    .rows()
                    .checked_mul(self.cols())
                    .ok_or_else(|| crate::Error::Invalid("matrix dimensions overflow".into()))?
        {
            return Err(crate::Error::Invalid(
                "matrix dimensions do not match buffers".into(),
            ));
        }
        Ok(())
    }
}
