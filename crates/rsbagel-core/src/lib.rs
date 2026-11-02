//! f64 in-memory BAGEL2 compatibility engine. See docs/network_profiles.md.
//! Translated behavior: Copyright (c) 2020 Hart Lab; MIT notice in docs/BAGEL2-LICENSE.txt.
pub mod error;
pub mod model;
pub mod profiling;
pub mod workers;
pub use error::{Error, Result};
pub use model::Matrix;
pub mod bf;
pub mod bf_types;
pub mod fold_change;
pub mod input;
pub mod kde;
pub mod resampling;
pub use bf_types::*;
pub use fold_change::{fc_to_bf_input, run_fc, FcConfig, FcResult};
pub mod multitarget;
pub mod network;
pub mod output;
pub mod precision_recall;
pub use precision_recall::{bf_to_pr_input, run_pr, PrRow};
pub fn run_bf(
    m: &Matrix,
    ess: &[String],
    non: &[String],
    cfg: &BfConfig,
    progress: Option<&(dyn Fn(usize, usize) -> std::ops::ControlFlow<()> + Sync)>,
) -> Result<BfResult> {
    let has_network = cfg.network.is_some() && !cfg.guide_level;
    let base_progress = |done, total| {
        progress.map_or(std::ops::ControlFlow::Continue(()), |f| {
            f(done, if has_network { total * 2 } else { total })
        })
    };
    let mut result = bf::run_bf_base(m, ess, non, cfg, Some(&base_progress))?;
    multitarget::apply_multitarget(&mut result, cfg)?;
    let main_count = result.partitions.len();
    let network_progress = |done, total| {
        progress.map_or(std::ops::ControlFlow::Continue(()), |f| {
            f(main_count + done, main_count + total)
        })
    };
    network::apply_network_with_progress(&mut result, ess, non, cfg, Some(&network_progress))?;
    Ok(result)
}
pub const SOURCE_COMMIT: &str = env!("RSBAGEL_SOURCE_COMMIT");
pub const SOURCE_DIRTY: &str = env!("RSBAGEL_SOURCE_DIRTY");
