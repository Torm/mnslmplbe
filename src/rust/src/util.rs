

/// Computes the log of the sum of exp applied to all numbers (the numbers should
/// be logarithmic values).
///
/// Uses the log sum exp trick to maintain the numerical stability.
pub fn log_sum_exp(logs: &[f64]) -> f64 {
    let mut max = f64::MIN;
    for l in logs {
        if *l > max {
            max = *l;
        }
    }
    let mut sum = 0.0;
    for l in logs {
        sum += (*l - max).exp();
    }
    max + sum.ln()
}

