//! Legacy NumPy RandomState plans and independently seeded native iteration plans.
use crate::bf_types::{BfConfig, Partition, Resampling, RngMode};
use crate::{Error, Result};
use rand::{RngCore, SeedableRng};
use rand_chacha::ChaCha8Rng;
use std::collections::HashSet;

#[derive(Clone)]
pub struct Mt19937 {
    state: [u32; 624],
    position: usize,
}
impl Mt19937 {
    pub fn new(seed: u32) -> Self {
        let mut state = [0; 624];
        state[0] = seed;
        for i in 1..624 {
            state[i] = 1812433253_u32
                .wrapping_mul(state[i - 1] ^ (state[i - 1] >> 30))
                .wrapping_add(i as u32);
        }
        Self {
            state,
            position: 624,
        }
    }
    pub fn next_u32(&mut self) -> u32 {
        if self.position == 624 {
            for i in 0..624 {
                let y = (self.state[i] & 0x80000000) | (self.state[(i + 1) % 624] & 0x7fffffff);
                self.state[i] = self.state[(i + 397) % 624]
                    ^ (y >> 1)
                    ^ if y & 1 != 0 { 0x9908b0df } else { 0 };
            }
            self.position = 0;
        }
        let mut y = self.state[self.position];
        self.position += 1;
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c5680;
        y ^= (y << 15) & 0xefc60000;
        y ^= y >> 18;
        y
    }
    pub fn uniform(&mut self) -> f64 {
        let a = self.next_u32() >> 5;
        let b = self.next_u32() >> 6;
        ((a as f64) * 67108864.0 + b as f64) / 9007199254740992.0
    }
    pub fn interval(&mut self, max: usize) -> usize {
        if max == 0 {
            return 0;
        }
        if max <= u32::MAX as usize {
            let mut mask = max as u32;
            mask |= mask >> 1;
            mask |= mask >> 2;
            mask |= mask >> 4;
            mask |= mask >> 8;
            mask |= mask >> 16;
            loop {
                let value = self.next_u32() & mask;
                if value as usize <= max {
                    return value as usize;
                }
            }
        }
        let mut mask = max as u64;
        mask |= mask >> 1;
        mask |= mask >> 2;
        mask |= mask >> 4;
        mask |= mask >> 8;
        mask |= mask >> 16;
        mask |= mask >> 32;
        loop {
            let value = (((self.next_u32() as u64) << 32) | self.next_u32() as u64) & mask;
            if value <= max as u64 {
                return value as usize;
            }
        }
    }
}

#[derive(Clone)]
pub struct Sampler {
    mode: RngMode,
    seed: u64,
    domain: u64,
    numpy: Mt19937,
    native: ChaCha8Rng,
}
fn splitmix(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9e3779b97f4a7c15);
    x = (x ^ (x >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94d049bb133111eb);
    x ^ (x >> 31)
}
impl Sampler {
    pub fn new(mode: RngMode, seed: u64) -> Result<Self> {
        if mode == RngMode::Numpy && seed > u32::MAX as u64 {
            return Err(Error::Invalid(
                "Seed must be between 0 and 2**32 - 1".into(),
            ));
        }
        Ok(Self {
            mode,
            seed,
            domain: 0,
            numpy: Mt19937::new(seed as u32),
            native: ChaCha8Rng::seed_from_u64(seed),
        })
    }
    fn begin_iteration(&mut self, iteration: usize, domain: u64) {
        if self.mode == RngMode::Rust {
            self.native = ChaCha8Rng::seed_from_u64(splitmix(
                self.seed
                    ^ splitmix(domain.wrapping_add(0xd1b54a32d192ed03))
                    ^ splitmix(iteration as u64),
            ));
        }
    }
    fn interval(&mut self, max: usize) -> usize {
        match self.mode {
            RngMode::Numpy => self.numpy.interval(max),
            RngMode::Rust => {
                if max == 0 {
                    return 0;
                }
                let bound = max as u64 + 1;
                let threshold = bound.wrapping_neg() % bound;
                loop {
                    let v = self.native.next_u64();
                    if v >= threshold {
                        return (v % bound) as usize;
                    }
                }
            }
        }
    }
    fn uniform(&mut self) -> f64 {
        match self.mode {
            RngMode::Numpy => self.numpy.uniform(),
            RngMode::Rust => (self.native.next_u64() >> 11) as f64 / 9007199254740992.0,
        }
    }
    pub fn training(
        &mut self,
        universe: &[usize],
        folds: usize,
        iterations: usize,
        bootstrap: bool,
    ) -> Result<Vec<Partition>> {
        let mut training = Training::new(universe, folds)?;
        (0..iterations)
            .map(|i| training.next_partition(self, bootstrap, i))
            .collect()
    }
}

// Rank selection preserves NumPy's changing interval bounds and the stable
// original-position order of Vec::remove without shifting every later entry.
// Small universes retain the Vec path; this is an implementation threshold only.
const NUMPY_RANK_BUCKET_MIN: usize = 4096;

struct RankBucket {
    tree: Vec<usize>,
    remaining: usize,
}
impl RankBucket {
    fn new(universe_len: usize, bucket: &[usize]) -> Self {
        let mut tree = vec![0; universe_len + 1];
        for &position in bucket {
            tree[position + 1] = 1;
        }
        for index in 1..tree.len() {
            let parent = index + index.isolate_lowest_one();
            if parent < tree.len() {
                let count = tree[index];
                tree[parent] += count;
            }
        }
        Self {
            tree,
            remaining: bucket.len(),
        }
    }
    fn remove_rank(&mut self, mut rank: usize) -> usize {
        debug_assert!(rank < self.remaining);
        let mut index = 0;
        let mut stride = 1;
        while stride <= (self.tree.len() - 1) / 2 {
            stride <<= 1;
        }
        while stride != 0 {
            let next = index + stride;
            if next < self.tree.len() && self.tree[next] <= rank {
                index = next;
                rank -= self.tree[next];
            }
            stride >>= 1;
        }
        let mut ancestor = index + 1;
        while ancestor < self.tree.len() {
            self.tree[ancestor] -= 1;
            ancestor += ancestor.isolate_lowest_one();
        }
        self.remaining -= 1;
        index
    }
    fn into_bucket(mut self) -> Vec<usize> {
        // Reverse the linear Fenwick construction. This rare conversion keeps
        // public Training calls correct if a caller switches to native mode.
        for index in (1..self.tree.len()).rev() {
            let parent = index + index.isolate_lowest_one();
            if parent < self.tree.len() {
                let count = self.tree[index];
                self.tree[parent] -= count;
            }
        }
        self.tree
            .into_iter()
            .skip(1)
            .enumerate()
            .filter_map(|(position, live)| (live != 0).then_some(position))
            .collect()
    }
}

pub struct Training {
    universe: Vec<usize>,
    bucket: Vec<usize>,
    rank_bucket: Option<RankBucket>,
    folds: usize,
    step: usize,
    domain: Option<u64>,
}
impl Training {
    pub fn new(universe: &[usize], folds: usize) -> Result<Self> {
        if folds == 0 {
            return Err(Error::Invalid(
                "cross-validation count must not be zero".into(),
            ));
        }
        Ok(Self {
            universe: universe.to_vec(),
            bucket: (0..universe.len()).collect(),
            rank_bucket: None,
            folds,
            step: 0,
            domain: None,
        })
    }
    pub fn next_partition(
        &mut self,
        sampler: &mut Sampler,
        bootstrap: bool,
        iteration: usize,
    ) -> Result<Partition> {
        let domain = *self.domain.get_or_insert_with(|| {
            let d = sampler.domain;
            sampler.domain = sampler.domain.wrapping_add(1);
            d
        });
        sampler.begin_iteration(iteration, domain);
        let n = self.universe.len();
        let mut selected = vec![false; n];
        if bootstrap {
            for _ in 0..n {
                let index = (sampler.uniform() * n as f64).floor() as usize;
                selected[index] = true;
            }
        } else {
            let bid = n / self.folds;
            if bid == 0 {
                return Err(Error::Legacy(format!(
                    "The number of genes is too small! n<{}",
                    self.folds
                )));
            }
            let draws = bid + usize::from(self.step < n % self.folds);
            if sampler.mode == RngMode::Numpy && n >= NUMPY_RANK_BUCKET_MIN {
                if self.rank_bucket.is_none() {
                    self.rank_bucket = Some(RankBucket::new(n, &self.bucket));
                    // No duplicate bucket storage in the ordinary NumPy path.
                    self.bucket = Vec::new();
                }
                let bucket = self.rank_bucket.as_mut().unwrap();
                for _ in 0..draws {
                    if bucket.remaining == 0 {
                        return Err(Error::Legacy("ValueError: high <= 0".into()));
                    }
                    let rank = sampler.interval(bucket.remaining - 1);
                    selected[bucket.remove_rank(rank)] = true;
                }
            } else {
                if let Some(bucket) = self.rank_bucket.take() {
                    self.bucket = bucket.into_bucket();
                }
                for _ in 0..draws {
                    if self.bucket.is_empty() {
                        return Err(Error::Legacy("ValueError: high <= 0".into()));
                    }
                    let index = sampler.interval(self.bucket.len() - 1);
                    selected[self.bucket.remove(index)] = true;
                }
            }
            self.step += 1;
        }
        let mut train = Vec::new();
        let mut test = Vec::new();
        for (i, &gene) in self.universe.iter().enumerate() {
            // Bootstrap selects training membership; CV selects withheld membership.
            if selected[i] == bootstrap {
                train.push(gene);
            } else {
                test.push(gene);
            }
        }
        Ok(Partition {
            iteration,
            train,
            test,
        })
    }
}

pub fn build_plans(
    genes: &[String],
    essentials: &[String],
    nonessentials: &[String],
    config: &BfConfig,
) -> Result<(Vec<Partition>, Sampler)> {
    let _profile_span = crate::profiling::span("sampling_plan");
    let mut sampler = Sampler::new(config.rng_mode, config.seed)?;
    let count = if config.no_resampling {
        1
    } else {
        match config.resampling {
            Resampling::CrossValidation => config.folds,
            Resampling::Bootstrap => config.bootstrap_iterations,
            Resampling::SmallSample => 100,
        }
    };
    if count == 0 {
        return Err(Error::Invalid(
            "iteration count must be greater than zero".into(),
        ));
    }
    if let Some(replay) = &config.replay {
        if replay.genes != genes {
            return Err(Error::Invalid(
                "replay gene identity/order does not match input".into(),
            ));
        }
        if replay.iterations.len() != count {
            return Err(Error::Invalid(format!(
                "replay iteration count {} does not match expected {count}",
                replay.iterations.len()
            )));
        }
        for (i, p) in replay.iterations.iter().enumerate() {
            validate_partition(p, genes.len(), i, config.no_resampling)?;
        }
        return Ok((replay.iterations.clone(), sampler));
    }
    let all: Vec<usize> = (0..genes.len()).collect();
    if config.no_resampling {
        return Ok((
            vec![Partition {
                iteration: 0,
                train: all.clone(),
                test: all,
            }],
            sampler,
        ));
    }
    let reference: HashSet<&str> = essentials
        .iter()
        .chain(nonessentials)
        .map(String::as_str)
        .collect();
    let universe: Vec<usize> = if config.resampling == Resampling::SmallSample {
        all.iter()
            .copied()
            .filter(|&i| reference.contains(genes[i].as_str()))
            .collect()
    } else {
        all.clone()
    };
    let mut plans = sampler.training(
        &universe,
        config.folds,
        count,
        config.resampling != Resampling::CrossValidation,
    )?;
    if config.resampling == Resampling::SmallSample {
        let others: Vec<usize> = all
            .into_iter()
            .filter(|&i| !reference.contains(genes[i].as_str()))
            .collect();
        for p in &mut plans {
            p.test.extend_from_slice(&others);
            p.test.sort_unstable();
            p.test.dedup();
        }
    }
    Ok((plans, sampler))
}

pub fn validate_partition(p: &Partition, n: usize, iteration: usize, overlap: bool) -> Result<()> {
    if p.iteration != iteration {
        return Err(Error::Invalid(
            "replay iteration IDs must be contiguous and ordered".into(),
        ));
    }
    let mut marks = vec![0u8; n];
    for (indices, bit) in [(&p.train, 1u8), (&p.test, 2u8)] {
        let mut previous = None;
        for &i in indices {
            if i >= n || previous.is_some_and(|v| v >= i) {
                return Err(Error::Invalid(
                    "replay indices must be unique, ordered and within gene range".into(),
                ));
            }
            marks[i] |= bit;
            previous = Some(i);
        }
    }
    if marks
        .iter()
        .any(|&m| if overlap { m != 3 } else { m != 1 && m != 2 })
    {
        return Err(Error::Invalid(
            "replay membership does not partition the input gene universe".into(),
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn sparse_rank_bucket_conversion_preserves_original_position_order() {
        let mut original = vec![1, 3, 4, 8, 12];
        let mut rank_bucket = RankBucket::new(16, &original);
        for rank in [0, 2, 1] {
            assert_eq!(rank_bucket.remove_rank(rank), original.remove(rank));
        }
        assert_eq!(rank_bucket.into_bucket(), original);
    }

    #[test]
    fn numpy_rank_bucket_matches_vec_memberships_and_rng_continuation() {
        // The unchanged Vec::remove algorithm is an independent oracle for rank
        // selection, including uneven folds and nonidentity universe ordering.
        for n in [4095, 4096, 4097, 8195] {
            let mut universe: Vec<usize> = (0..n).rev().map(|i| i * 3 + 7).collect();
            universe.rotate_left(137);
            for folds in [3, 10, 17] {
                for seed in [0, 42, u32::MAX as u64] {
                    let mut sampler = Sampler::new(RngMode::Numpy, seed).unwrap();
                    let mut oracle = Sampler::new(RngMode::Numpy, seed).unwrap();
                    let mut training = Training::new(&universe, folds).unwrap();
                    let mut bucket: Vec<usize> = (0..n).collect();
                    for iteration in 0..folds {
                        let draws = n / folds + usize::from(iteration < n % folds);
                        let mut selected = vec![false; n];
                        for _ in 0..draws {
                            let rank = oracle.interval(bucket.len() - 1);
                            selected[bucket.remove(rank)] = true;
                        }
                        let expected_train: Vec<usize> = universe
                            .iter()
                            .enumerate()
                            .filter_map(|(i, &gene)| (!selected[i]).then_some(gene))
                            .collect();
                        let expected_test: Vec<usize> = universe
                            .iter()
                            .enumerate()
                            .filter_map(|(i, &gene)| selected[i].then_some(gene))
                            .collect();
                        let actual = training
                            .next_partition(&mut sampler, false, iteration)
                            .unwrap();
                        assert_eq!(actual.iteration, iteration);
                        assert_eq!(
                            actual.train, expected_train,
                            "n={n}, folds={folds}, seed={seed}"
                        );
                        assert_eq!(
                            actual.test, expected_test,
                            "n={n}, folds={folds}, seed={seed}"
                        );
                        assert_eq!(training.rank_bucket.is_some(), n >= NUMPY_RANK_BUCKET_MIN);
                        assert_eq!(sampler.numpy.position, oracle.numpy.position);
                        assert_eq!(sampler.numpy.state, oracle.numpy.state);
                    }
                    assert!(bucket.is_empty());
                    assert!(training.next_partition(&mut sampler, false, folds).is_err());
                    // Rejection-sampled integers and doubles must continue at the
                    // identical state, rather than merely produce the same fold.
                    for bound in [0, 1, 17, 4096, u32::MAX as usize] {
                        assert_eq!(sampler.interval(bound), oracle.interval(bound));
                    }
                    for _ in 0..8 {
                        assert_eq!(sampler.uniform().to_bits(), oracle.uniform().to_bits());
                    }
                }
            }
        }
    }

    #[test]
    fn numpy_words_and_doubles() {
        let mut m = Mt19937::new(42);
        assert_eq!(m.next_u32(), 1608637542);
        assert_eq!(m.next_u32(), 3421126067);
        let mut m = Mt19937::new(42);
        assert_eq!(m.uniform(), 0.3745401188473625);
        assert_eq!(m.uniform(), 0.9507143064099162);
        let mut m = Mt19937::new(42);
        assert_eq!(m.interval(0), 0);
        assert_eq!(m.next_u32(), 1608637542);
    }
    #[test]
    fn numpy_cv_fixture() {
        let mut s = Sampler::new(RngMode::Numpy, 42).unwrap();
        let ps = s
            .training(&(0..12).collect::<Vec<_>>(), 3, 3, false)
            .unwrap();
        assert_eq!(ps[0].test, vec![3, 5, 6, 9]);
        assert_eq!(ps.iter().map(|p| p.test.len()).sum::<usize>(), 12);
    }
    #[test]
    #[ignore = "HAP1 regression fixture is not distributed in the snapshot"]
    fn hap1_numpy_plan_matches_reference() {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/golden/hap1-cv-final/replay.json");
        let replay: crate::bf_types::ReplayManifest =
            serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap();
        // The same defaults are used for direct Rust callers and omitted C-ABI JSON fields.
        let config: BfConfig = serde_json::from_str("{}").unwrap();
        assert_eq!(config.rng_mode, RngMode::Numpy);
        assert_eq!(RngMode::default(), RngMode::Numpy);
        assert_eq!(config.network_profile, crate::NetworkProfile::Legacy);
        let explicit_native: BfConfig = serde_json::from_str(r#"{"rng_mode":"rust"}"#).unwrap();
        assert_eq!(explicit_native.rng_mode, RngMode::Rust);
        let (actual, _) = build_plans(&replay.genes, &[], &[], &config).unwrap();
        assert_eq!(
            serde_json::to_vec(&actual).unwrap(),
            serde_json::to_vec(&replay.iterations).unwrap()
        );
    }
    #[test]
    #[ignore = "explicit native G7 resampling-plan diagnosis"]
    fn capture_bootstrap_plan_hashes() {
        use sha2::{Digest, Sha256};
        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
        let replay: crate::bf_types::ReplayManifest = serde_json::from_str(
            &std::fs::read_to_string(root.join("tests/golden/hap1-cv-final/replay.json")).unwrap(),
        )
        .unwrap();
        let output = root.join("tests/external/qualification/native-rng");
        for (mode, name) in [
            (RngMode::Numpy, "numpy-compatible"),
            (RngMode::Rust, "native"),
        ] {
            let config = BfConfig {
                rng_mode: mode,
                resampling: Resampling::Bootstrap,
                bootstrap_iterations: 1000,
                ..Default::default()
            };
            let (plans, _) = build_plans(&replay.genes, &[], &[], &config).unwrap();
            let hashes:Vec<serde_json::Value>=plans.iter().map(|p| {
                let mut hash=Sha256::new();
                for i in [p.iteration,p.train.len(),p.test.len()].into_iter().chain(p.train.iter().copied()).chain(p.test.iter().copied()) {hash.update((i as u64).to_le_bytes());}
                serde_json::json!({"iteration":p.iteration,"train_count":p.train.len(),"test_count":p.test.len(),"sha256":format!("{:x}",hash.finalize())})
            }).collect();
            let file = std::fs::OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(output.join(format!("{name}-partition-hashes.json")))
                .unwrap();
            serde_json::to_writer_pretty(file, &hashes).unwrap();
        }
    }
    #[test]
    fn cv_exhaustion_is_retained() {
        let mut s = Sampler::new(RngMode::Numpy, 42).unwrap();
        let mut t = Training::new(&(0..120).collect::<Vec<_>>(), 10).unwrap();
        for i in 0..10 {
            assert_eq!(t.next_partition(&mut s, false, i).unwrap().test.len(), 12);
        }
        assert!(t.next_partition(&mut s, false, 10).is_err());
    }
}
