//! Opt-in exploratory instrumentation. No arithmetic or scientific state is changed.
//! Durations are inclusive wall spans; overlapping worker spans must not be summed
//! as elapsed time. Enable only with RSBAGEL_PROFILE pointing to a JSONL sink.
use std::collections::BTreeMap;
use std::io::Write;
use std::sync::{Mutex, OnceLock};
use std::time::Instant;

#[derive(Default, serde::Serialize)]
struct Total {
    calls: u64,
    inclusive_wall_seconds: f64,
}
struct State {
    path: Option<std::path::PathBuf>,
    totals: Mutex<BTreeMap<&'static str, Total>>,
}
fn state() -> &'static State {
    static STATE: OnceLock<State> = OnceLock::new();
    STATE.get_or_init(|| State {
        path: std::env::var_os("RSBAGEL_PROFILE").map(Into::into),
        totals: Mutex::new(BTreeMap::new()),
    })
}
pub struct Span {
    stage: &'static str,
    started: Option<Instant>,
}
pub fn span(stage: &'static str) -> Span {
    Span {
        stage,
        started: state().path.as_ref().map(|_| Instant::now()),
    }
}
impl Drop for Span {
    fn drop(&mut self) {
        if let Some(start) = self.started {
            let seconds = start.elapsed().as_secs_f64();
            if let Ok(mut totals) = state().totals.lock() {
                let total = totals.entry(self.stage).or_default();
                total.calls += 1;
                total.inclusive_wall_seconds += seconds;
            }
        }
    }
}
/// Writes a cumulative snapshot at scope exit. Profiling I/O failures never alter
/// scientific results; the supervisor requires the sink and rejects missing data.
pub struct Session;
impl Default for Session {
    fn default() -> Self {
        Self::new()
    }
}
impl Session {
    pub fn new() -> Self {
        let _ = state();
        Self
    }
}
impl Drop for Session {
    fn drop(&mut self) {
        if let Some(path) = &state().path {
            if let Ok(totals) = state().totals.lock() {
                if let Ok(mut file) = std::fs::OpenOptions::new()
                    .create(true)
                    .append(true)
                    .open(path)
                {
                    let record = serde_json::json!({"schema":"rsbagel-profile-spans-v1","pid":std::process::id(),"source_commit":crate::SOURCE_COMMIT,"source_dirty":crate::SOURCE_DIRTY,"stages":*totals,"semantics":"cumulative inclusive wall spans; worker sums overlap; CPU and RSS are supervisor process metrics"});
                    let _ = serde_json::to_writer(&mut file, &record);
                    let _ = writeln!(file);
                }
            }
        }
    }
}
