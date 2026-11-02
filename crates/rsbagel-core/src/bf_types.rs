use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "kebab-case")]
pub enum Resampling {
    #[default]
    CrossValidation,
    Bootstrap,
    SmallSample,
}
#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "lowercase")]
pub enum RngMode {
    /// Deterministic native ChaCha8; experimental and not production-qualified.
    Rust,
    /// BAGEL2-compatible NumPy RandomState sampling; production default.
    #[default]
    Numpy,
}
#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "lowercase")]
pub enum NetworkProfile {
    #[default]
    Legacy,
    Corrected,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Partition {
    pub iteration: usize,
    pub train: Vec<usize>,
    pub test: Vec<usize>,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ReplayManifest {
    pub genes: Vec<String>,
    pub iterations: Vec<Partition>,
    #[serde(default)]
    pub network_iterations: Vec<Partition>,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct AlignmentRecord {
    pub guide: String,
    pub perfect_loci: usize,
    pub perfect_genes: usize,
    pub mismatch_loci: usize,
    pub mismatch_genes: usize,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(default)]
pub struct BfConfig {
    /// Zero-based numeric matrix columns. Repeats intentionally retained for scoring.
    pub columns: Vec<usize>,
    pub column_labels: Vec<String>,
    pub resampling: Resampling,
    pub no_resampling: bool,
    pub folds: usize,
    pub bootstrap_iterations: usize,
    pub seed: u64,
    pub rng_mode: RngMode,
    pub threads: usize,
    pub replay: Option<ReplayManifest>,
    pub guide_level: bool,
    pub equalise_guides: Option<f64>,
    pub equalise_reps: Option<f64>,
    pub alignment: Option<Vec<AlignmentRecord>>,
    pub max_perfect_loci: usize,
    pub max_mismatch_loci: usize,
    pub network: Option<Vec<(String, String)>>,
    pub network_profile: NetworkProfile,
    pub dump_stages: bool,
}
impl Default for BfConfig {
    fn default() -> Self {
        Self {
            columns: vec![0],
            column_labels: vec![],
            resampling: Resampling::CrossValidation,
            no_resampling: false,
            folds: 10,
            bootstrap_iterations: 1000,
            seed: 42,
            rng_mode: RngMode::Numpy,
            threads: crate::workers::auto_workers(),
            replay: None,
            guide_level: false,
            equalise_guides: None,
            equalise_reps: None,
            alignment: None,
            max_perfect_loci: 10,
            max_mismatch_loci: 10,
            network: None,
            network_profile: NetworkProfile::Legacy,
            dump_stages: false,
        }
    }
}
#[derive(Clone, Debug, Serialize)]
pub struct GeneScore {
    pub gene: String,
    pub raw_mean: f64,
    pub mean: f64,
    pub std: f64,
    pub count: usize,
    pub normalized: Option<f64>,
    pub boosted: Option<f64>,
    pub boosted_std: Option<f64>,
    pub boosted_normalized: Option<f64>,
    pub network_evidence: Option<bool>,
    pub network_contribution_count: Option<usize>,
    pub network_zero_applied: Option<bool>,
}
#[derive(Clone, Debug, Serialize)]
pub struct GuideScore {
    pub guide: String,
    pub gene: String,
    pub mean_replicates: Vec<f64>,
    pub std_replicates: Vec<f64>,
    pub total: f64,
    pub count: usize,
}
#[derive(Clone, Debug, Serialize)]
pub struct IterationTrace {
    pub partition: Partition,
    pub essential: Vec<f64>,
    pub nonessential: Vec<f64>,
    pub essential_factor: f64,
    pub nonessential_factor: f64,
    pub essential_covariance: f64,
    pub nonessential_covariance: f64,
    pub grid: Vec<f64>,
    pub nonessential_density: Vec<f64>,
    pub support_indices: Vec<usize>,
    pub subgrid: Vec<f64>,
    pub logratios: Vec<f64>,
    pub minimum_indices: Vec<usize>,
    pub xmin: f64,
    pub xmax: f64,
    pub regression_x: Vec<f64>,
    pub regression_y: Vec<f64>,
    pub slope: f64,
    pub intercept: f64,
}
#[derive(Clone, Debug, Serialize)]
pub struct BfResult {
    /// First-seen gene order, independent from lexicographic output serialization.
    pub genes: Vec<GeneScore>,
    pub guides: Vec<GuideScore>,
    pub gene_guides: Vec<Vec<usize>>,
    pub partitions: Vec<Partition>,
    /// Guide, iteration participation, replicate. Retained for raw qualification.
    pub guide_iterations: Vec<Vec<Vec<f64>>>,
    pub traces: Vec<IterationTrace>,
    pub network_partitions: Vec<Partition>,
    pub network_scores: Vec<Option<f64>>,
    pub network_events: Vec<NetworkEvent>,
    pub network_traces: Vec<NetworkTrace>,
    pub diagnostics: Vec<String>,
    pub correction_coefficients: Option<[f64; 2]>,
    pub correction_trace: Option<MultiTargetTrace>,
}

#[derive(Clone, Debug, Serialize)]
pub struct NetworkEvent {
    pub gene: String,
    pub score: f64,
    pub contribution: f64,
}
#[derive(Clone, Debug, Serialize)]
pub struct NetworkTrace {
    pub partition: Partition,
    pub essential: Vec<f64>,
    pub nonessential: Vec<f64>,
    pub xmin: f64,
    pub xmax: f64,
    pub regression_x: Vec<f64>,
    pub regression_y: Vec<f64>,
    pub slope: f64,
    pub intercept: f64,
}

#[derive(Clone, Debug, Serialize)]
pub struct MultiTargetTrace {
    pub guides: Vec<String>,
    pub x: Vec<[f64; 2]>,
    pub y: Vec<f64>,
    pub coefficients_unclipped: [f64; 2],
    pub coefficients_clipped: [f64; 2],
    pub intercept: f64,
}
