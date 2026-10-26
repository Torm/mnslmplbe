use std::cmp::{min, Ordering};
use std::ops::{Deref, DerefMut};
use std::thread::available_parallelism;
use console::Term;
use dashmap::DashMap;
use ndarray::Array2;
use parking_lot::Mutex;
use rayon::slice::ParallelSliceMut;
use statrs::function::factorial::binomial;
use statrs::function::gamma::ln_gamma;
use crate::iter::CombinationIter;
use crate::util::log_sum_exp;


/// Blanket selection configuration.
#[derive(Copy, Clone)]
pub enum BlanketSelect { All, Cap(f64), CombinationLimit(f64) }


/// Blanket expansion method with progress terminal output. // TODO Support RStudio terminal better
pub fn be_progress(max_k: u32, blanket_select_from_k: u32, domain: &[u32], sample: &Array2<u32>, alpha: f64, criteria: BlanketSelect, terminal: Option<&Term>) -> Result<Array2<f64>, String> {
    let n = domain.len();

    let mut blanket_member_posteriors = Array2::from_shape_simple_fn((n, n), || 0.0);

    if let Some(terminal) = terminal {
        terminal.write_line(&format!("MN Network Α Posteriors k={} α={} n={} m={}", max_k, alpha, domain.len(), sample.nrows())).unwrap();
        terminal.write_line(&format!("...")).unwrap();
    }

    for i in 0..n {
        if let Some(terminal) = terminal {
            terminal.move_cursor_left(10000).unwrap();
            terminal.move_cursor_up(1).unwrap();
            terminal.write_line(&format!("Node {}/{}...", i + 1, n)).unwrap();
        }
        let posterior = be_node_posteriors(sample, max_k, blanket_select_from_k, i as u32, domain, alpha, criteria, terminal).unwrap();
        for (p, s) in posterior {
            blanket_member_posteriors[(i, p as usize)] = s;
        }
    }
    if let Some(terminal) = terminal {
        terminal.move_cursor_left(10000).unwrap();
        terminal.move_cursor_up(1).unwrap();
    }

    Ok(blanket_member_posteriors)
}


/// Compute the blanket member posteriors for a node.
pub fn be_node_posteriors(
    sample: &Array2<u32>, max_k: u32, blanket_select_from_k: u32, u: u32, domain: &[u32], alpha: f64, criteria: BlanketSelect, terminal: Option<&Term>
) -> Result<Box<[(u32, f64)]>, String> {
    let n = domain.len();

    let mut selected_blanket = vec![]; // Nodes investigated. Initially all nodes are considered possible blanket members.
    let mut sorted_posteriors = Vec::with_capacity(n as usize);
    for i in 0..n {
        sorted_posteriors.push((i as u32, 0.0f64));
    }
    let mut sorted_posteriors = sorted_posteriors.into_boxed_slice();

    let mut log_blanket_scores = vec![f64::MIN; n as usize]; // Accumulated scores for each blanket. Shared across iterations. Numerator each in posterior calculation.
    let mut log_total_blanket_scores = f64::MIN; // Total accumulated blanket score. Denominator in each posterior calculation.

    let mut scores_mutex = Mutex::new((log_blanket_scores, log_total_blanket_scores));

    if let Some(terminal) = terminal {
        terminal.write_line(&format!("...")).unwrap();
    }


    for k in 1..=max_k {

        // Select blankets based on criteria.
        selected_blanket.clear();
        // Do not start blanket filtering until k = blanket_select_from_k.
        let criteria = if k < blanket_select_from_k {
            BlanketSelect::All
        } else {
            criteria
        };
        match criteria {
            BlanketSelect::All => { // No need to do anything; all nodes remain in the blanket.
                for (p, s) in sorted_posteriors.iter() {
                    if *p == u {
                        continue;
                    }
                    selected_blanket.push(*p);
                }
            }
            BlanketSelect::Cap(factor) => {
                let min_posterior = sorted_posteriors[0].1 * factor;
                for (p, s) in sorted_posteriors.iter() {
                    if *p == u {
                        continue;
                    }
                    if *s >= min_posterior {
                        selected_blanket.push(*p);
                    } else {
                        break;
                    }
                }
            }
            BlanketSelect::CombinationLimit(s) => {
                let np = min(largest_combination_lesser_than(k as u64, s), sorted_posteriors.len() as u64);
                for i in 0..np {
                    let (p, _) = sorted_posteriors[i as usize];
                    if p == u {
                        continue;
                    }
                    selected_blanket.push(p);
                }
            }
        }

        if let Some(terminal) = terminal {
            terminal.move_cursor_left(10000).unwrap();
            terminal.move_cursor_up(1).unwrap();
            //terminal.write_line(&format!("Scoring blankets of size {}, selected {} members, {} combinations...", k, selected_blankets.len(), binomial(selected_blankets.len() as u64, k as u64))).unwrap();/// ///TODO Enable
        }

        // Iterate over the blankets that are possible given the selected variables.
        let mut combination_iter = Mutex::new(CombinationIter::new(selected_blanket.len(), k as usize));
        let threads = available_parallelism().unwrap().into();
        std::thread::scope(|s| {
            for _ in 0..threads {
                s.spawn(|| {
                    // Each thread iterates over the next available blanket combination.
                    let mut blanket_buffer = vec![0; k as usize].into_boxed_slice(); // Holds the blanket k-combination. Initially invalid.
                    loop {
                        // Iterate to the next k-combination of the selected blankets.
                        let mut locked_iter = combination_iter.lock();
                        if !locked_iter.next() { // TODO: Check num of selected members >= k
                            break;
                        }
                        // Read the blanket combination into the buffer.
                        for i in 0..blanket_buffer.len() {
                            blanket_buffer[i] = selected_blanket[locked_iter.at()[i] as usize];
                        }
                        drop(locked_iter);
                        // Compute the log blanket score.
                        let lfs = log_blanket_score(sample, u, &blanket_buffer, domain, alpha).unwrap();
                        // Add the log-blanket-scores to relevant neighbours and total. Use log-sum-exp to maintain numerical stability.
                        let mut locked_scores = scores_mutex.lock();
                        for p in blanket_buffer.iter() {
                            let current_b_score = locked_scores.deref().0[*p as usize];
                            locked_scores.deref_mut().0[*p as usize] = log_sum_exp(&[current_b_score, lfs]);
                        }
                        let current_total = locked_scores.deref().1;
                        locked_scores.deref_mut().1 = log_sum_exp(&[current_total, lfs]);
                        drop(locked_scores);
                    }
                });
            }
        });

        {
            let (log_b_scores, log_total_b_scores) = scores_mutex.get_mut();

            for i in 0..n {
                sorted_posteriors[i] = (i as u32, log_b_scores[i] - *log_total_b_scores);
            }
        }

        sorted_posteriors.par_sort_unstable_by(|&(_, s1), &(_, s2)| {
            if let Some(c) = s2.partial_cmp(&s1) {
                c
            } else {
                eprintln!("None partial cmp of {} and {}", s1, s2);
                Ordering::Equal
            }

        }); // Sort the posteriors.

    }

    if let Some(terminal) = terminal {
        terminal.move_cursor_left(10000).unwrap();
        terminal.move_cursor_up(1).unwrap();
    }

    Ok(sorted_posteriors)
}


/// Compute the log blanket score of a blanket.
///
/// Counts the configurations of the blanket by iterating through the sample.
pub fn log_blanket_score(sample: &Array2<u32>, node: u32, blanket: &[u32], domain: &[u32], α: f64) -> Result<f64, String> {
    let bn = blanket.len();
    let m = sample.nrows();
    let n = domain.len();
    // Count the number of joint and blanket-only configs in the data sample.
    let counts: DashMap<Box<[u32]>, u32> = DashMap::new();
    let mut config_buffer = vec![0; bn + 1].into_boxed_slice(); // Buffer for the blanket config, and the joint's node config in the last slot.
    for i in 0..m {
        // Read the config in the current row into the buffer.
        for j in 0..bn {
            config_buffer[j] = sample[(i, blanket[j] as usize)];
        }
        config_buffer[bn] = sample[(i, node as usize)];
        // Increment the joint config count.
        if let Some(mut count) = counts.get_mut(&config_buffer) {
            *count += 1;
        } else {
            counts.insert(config_buffer.clone(), 1);
        }
        // Increment the blanket config count.
        let blanket_config_sub_buffer = &config_buffer[0..blanket.len()];
        if let Some(mut count) = counts.get_mut(blanket_config_sub_buffer) {
            *count += 1;
        } else {
            counts.insert(Box::from(blanket_config_sub_buffer), 1);
        }
    }
    // Compute the total number of configurations of the blanket.
    let mut blanket_configs = 1;
    for p in &config_buffer[0..bn] { // Compute the number of blanket configs.
        blanket_configs *= domain[*p as usize];
    }
    // Compute blanket and joint α. The BDe prior distributes α evenly over the blanket configurations. // TODO: More prior configurations than BDe
    let blanket_term_α = α / blanket_configs as f64;
    let joint_term_α = blanket_term_α / domain[node as usize] as f64;
    // Compute the log blanket score for the blanket and subject node.
    let mut log_blanket_score = 0.0;
    for entry in counts.iter() {
        let config = entry.key();
        let count = *entry.value();
        if config.len() == bn + 1 { // Joint config: contains both the blanket and the variable.
            log_blanket_score += log_blanket_score_joint_config_term(count, joint_term_α);
        } else if config.len() == bn { // Blanket config: contains only values for the blanket.
            log_blanket_score += log_blanket_score_blanket_config_term(count, blanket_term_α);
        } else {
            unreachable!();
        }
    }
    Ok(log_blanket_score)
}


pub fn log_blanket_score_blanket_config_term(blanket_count: u32, blanket_alpha: f64) -> f64 {
    ln_gamma(blanket_alpha) - ln_gamma(blanket_alpha + blanket_count as f64)
}


pub fn log_blanket_score_joint_config_term(joint_count: u32, joint_alpha: f64) -> f64 {
    ln_gamma(joint_alpha + joint_count as f64) - ln_gamma(joint_alpha)
}


fn largest_combination_lesser_than(k: u64, s: f64) -> u64 {
    let mut n = k;
    loop {
        n += 1;
        let c = binomial(n, k);
        if c > s {
            return n - 1;
        }
    }
}
