fn main() {
    let output = std::process::Command::new("git")
        .args(["rev-parse", "HEAD"])
        .output();
    let commit = output
        .ok()
        .filter(|o| o.status.success())
        .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
        .unwrap_or_else(|| "source-unavailable".into());
    println!("cargo:rustc-env=RSBAGEL_SOURCE_COMMIT={commit}");
    let dirty = std::process::Command::new("git")
        .args(["status", "--porcelain", "--untracked-files=no"])
        .output()
        .map(|o| !o.stdout.is_empty())
        .unwrap_or(true);
    println!("cargo:rustc-env=RSBAGEL_SOURCE_DIRTY={dirty}");
    println!("cargo:rerun-if-changed=../../.git/HEAD");
    println!("cargo:rerun-if-changed=../../.git/refs/heads/codex/rsbagel-campaign");
}
