//! Scalar f64 kernels matching the pinned NumPy/SciPy BAGEL2 path.
use crate::{Error, Result};

/// NumPy's pairwise contiguous f64 reducer (eight accumulators, 128 element blocks).
pub fn numpy_sum(x: &[f64]) -> f64 {
    let n = x.len();
    if n < 8 {
        x.iter().fold(-0.0, |a, b| a + b)
    } else if n <= 128 {
        let mut r = [0.0; 8];
        r.copy_from_slice(&x[..8]);
        let mut i = 8;
        while i + 8 <= n {
            for j in 0..8 {
                r[j] += x[i + j];
            }
            i += 8;
        }
        let mut sum = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
        for value in &x[i..] {
            sum += value;
        }
        sum
    } else {
        let half = (n / 2) & !7;
        numpy_sum(&x[..half]) + numpy_sum(&x[half..])
    }
}
pub fn mean(x: &[f64]) -> f64 {
    numpy_sum(x) / x.len() as f64
}
pub fn std(x: &[f64]) -> f64 {
    let m = mean(x);
    let deviations: Vec<f64> = x.iter().map(|v| (v - m) * (v - m)).collect();
    (numpy_sum(&deviations) / x.len() as f64).sqrt()
}
pub fn round_hundredth(x: f64) -> f64 {
    (x * 100.0).round_ties_even() / 100.0
}

/// NumPy uses ceil length and the representable `(start+step)-start` stride.
pub fn arange(start: f64, stop: f64, step: f64) -> Result<Vec<f64>> {
    if !start.is_finite() || !stop.is_finite() || !step.is_finite() || step == 0.0 {
        return Err(Error::Numerical(
            "arange: non-finite bounds or zero step".into(),
        ));
    }
    let length = ((stop - start) / step).ceil();
    if length <= 0.0 {
        return Ok(Vec::new());
    }
    if !length.is_finite() || length > isize::MAX as f64 / 8.0 {
        return Err(Error::Numerical(
            "arange: maximum allowed size exceeded".into(),
        ));
    }
    let length = length as usize;
    let mut out = Vec::new();
    out.try_reserve_exact(length)
        .map_err(|_| Error::Numerical("arange: allocation failed".into()))?;
    let delta = (start + step) - start;
    for i in 0..length {
        out.push(if i == 0 {
            start
        } else if i == 1 {
            start + step
        } else {
            start + i as f64 * delta
        });
    }
    Ok(out)
}

#[derive(Clone, Debug)]
pub struct Kde {
    pub dataset: Vec<f64>,
    pub factor: f64,
    pub covariance: f64,
    /// Cholesky scale, deliberately evaluated separately from sqrt(covariance).
    pub bandwidth: f64,
    whitened: Vec<f64>,
    weighted_normalizer: f64,
}
impl Kde {
    pub fn new(dataset: &[f64]) -> Result<Self> {
        let _profile_span = crate::profiling::span("kde_fit");
        if dataset.len() < 2 {
            return Err(Error::Numerical(
                "`dataset` input should have multiple elements.".into(),
            ));
        }
        if dataset.iter().any(|v| !v.is_finite()) {
            return Err(Error::Numerical(
                "array must not contain infs or NaNs".into(),
            ));
        }
        let n = dataset.len();
        let weight = 1.0 / n as f64;
        let weights = vec![weight; n];
        let weight_squares = vec![weight * weight; n];
        let wsum = numpy_sum(&weights);
        let w2sum = numpy_sum(&weight_squares);
        let neff = 1.0 / w2sum;
        let factor = neff.powf(-0.2);
        let weighted: Vec<f64> = dataset.iter().map(|x| x * weight).collect();
        let center = numpy_sum(&weighted) / wsum;
        // np.cov(aweights=weights) evaluates dot(centered, centered*weights).
        let mut covariance_sum = 0.0;
        for x in dataset {
            let d = x - center;
            covariance_sum += d * (d * weight);
        }
        let data_covariance = covariance_sum / (wsum - w2sum / wsum);
        if !data_covariance.is_finite() || data_covariance <= 0.0 {
            return Err(Error::Numerical(
                "singular data covariance matrix in gaussian_kde".into(),
            ));
        }
        let covariance = data_covariance * factor.powi(2);
        let bandwidth = data_covariance.sqrt() * factor;
        let whitened = dataset.iter().map(|x| x / bandwidth).collect();
        let weighted_normalizer = weight * ((2.0 * std::f64::consts::PI).powf(-0.5) / bandwidth);
        Ok(Self {
            dataset: dataset.to_vec(),
            factor,
            covariance,
            bandwidth,
            whitened,
            weighted_normalizer,
        })
    }
    pub fn evaluate(&self, point: f64) -> f64 {
        let point = point / self.bandwidth;
        let mut sum = 0.0;
        for observation in &self.whitened {
            let difference = observation - point;
            sum += (-0.5 * difference * difference).exp() * self.weighted_normalizer;
        }
        sum
    }
    pub fn evaluate_many(&self, points: &[f64]) -> Vec<f64> {
        let _profile_span = crate::profiling::span("kde_evaluate_many");
        points.iter().map(|&x| self.evaluate(x)).collect()
    }
}

/// Only slope and intercept are consumed by BAGEL2. Singleton input yields NaNs,
/// as pinned SciPy does; callers implementing corrected semantics check finiteness.
pub fn linregress(x: &[f64], y: &[f64]) -> Result<(f64, f64)> {
    if x.is_empty() || x.len() != y.len() {
        return Err(Error::Numerical(
            "Inputs must not be empty and must have equal lengths.".into(),
        ));
    }
    if x.len() > 1 && x.iter().all(|v| *v == x[0]) {
        return Err(Error::Numerical(
            "Cannot calculate a linear regression if all x values are identical".into(),
        ));
    }
    let xm = mean(x);
    let ym = mean(y);
    let mut xx = 0.0;
    let mut xy = 0.0;
    for (&a, &b) in x.iter().zip(y) {
        xx += (a - xm) * (a - xm);
        xy += (a - xm) * (b - ym);
    }
    let n = x.len() as f64;
    let slope = (xy / n) / (xx / n);
    Ok((slope, ym - slope * xm))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn numpy_grid_stride_and_rounding() {
        let x = arange(-10.0, 2.0, 0.01).unwrap();
        assert_eq!(x.len(), 1200);
        assert_eq!(x[1].to_bits(), (-9.99_f64).to_bits());
        assert_eq!(x[1199].to_bits(), (1.9899999999997444_f64).to_bits());
        assert_eq!(round_hundredth(0.125), 0.12);
        assert_eq!(round_hundredth(-0.125), -0.12);
    }
    #[test]
    fn degeneracies_are_not_repaired() {
        assert!(Kde::new(&[1.0]).is_err());
        assert!(Kde::new(&[1.0, 1.0]).is_err());
        assert!(Kde::new(&[0.0, f64::INFINITY]).is_err());
        assert!(linregress(&[2.0, 2.0], &[1.0, 3.0]).is_err());
        assert!(linregress(&[1.0], &[2.0]).unwrap().0.is_nan());
        assert!(mean(&[]).is_nan());
    }
    #[test]
    fn simple_kde_and_regression() {
        let k = Kde::new(&[-1.0, 0.0, 1.0]).unwrap();
        assert!((k.factor - 0.8027415617602307).abs() < 1e-15);
        assert!((k.evaluate(0.0) - 0.31815622092640145).abs() < 1e-14);
        let (a, b) = linregress(&[0.0, 1.0, 2.0], &[1.0, 3.0, 5.0]).unwrap();
        assert_eq!((a, b), (2.0, 1.0));
    }
}
