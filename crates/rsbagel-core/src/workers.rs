//! Bounded automatic worker selection. No workload-size heuristic.

/// Resolve a topology snapshot; zero/unavailable topology falls back to parallelism.
pub fn resolve_auto(performance_cores: Option<usize>, available: usize) -> usize {
    performance_cores
        .filter(|n| *n > 0)
        .unwrap_or(available)
        .clamp(1, 8)
}

pub fn auto_workers() -> usize {
    let available = std::thread::available_parallelism().map_or(1, usize::from);
    resolve_auto(performance_cores(), available)
}

#[cfg(all(target_os = "macos", target_arch = "aarch64"))]
fn performance_cores() -> Option<usize> {
    // Darwin perflevel0 is the performance cluster on supported Apple Silicon.
    let output = std::process::Command::new("/usr/sbin/sysctl")
        .args(["-n", "hw.perflevel0.physicalcpu"])
        .output()
        .ok()?;
    output.status.success().then_some(())?;
    std::str::from_utf8(&output.stdout)
        .ok()?
        .trim()
        .parse()
        .ok()
}

#[cfg(not(all(target_os = "macos", target_arch = "aarch64")))]
fn performance_cores() -> Option<usize> {
    None
}

/// Parse the CLI policy while preserving explicit worker counts without a cap.
pub fn parse_workers(value: &str) -> std::result::Result<usize, String> {
    if value == "auto" {
        return Ok(auto_workers());
    }
    value
        .parse::<usize>()
        .ok()
        .filter(|n| *n > 0)
        .ok_or_else(|| "workers must be auto or a positive integer".into())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn topology_and_fallback() {
        assert_eq!(resolve_auto(Some(6), 12), 6);
        assert_eq!(resolve_auto(Some(12), 16), 8);
        assert_eq!(resolve_auto(None, 16), 8);
        assert_eq!(resolve_auto(Some(0), 4), 4);
        assert_eq!(resolve_auto(None, 0), 1);
    }
    #[test]
    fn explicit_and_invalid() {
        assert_eq!(parse_workers("1").unwrap(), 1);
        assert_eq!(parse_workers("24").unwrap(), 24);
        assert!((1..=8).contains(&parse_workers("auto").unwrap()));
        for value in ["0", "-1", "1.5", "AUTO", "", "999999999999999999999999"] {
            assert!(parse_workers(value).is_err(), "{value}");
        }
    }
}
