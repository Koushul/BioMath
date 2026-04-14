//! Johnson et al. (*Cell* 2025) hypothesis grammar: CSV rows and canonical English sentences.
//!
//! Extended CSV schema (semicolon-separated, matching the methods supplement):
//! - **8 columns** (legacy): `cell_type ; signal ; response ; behavior ; max_response ; half_max ; hill_power ; applies_to_dead` — **Hill** response only.
//! - **11 columns**: append `response_kind ; param_a ; param_b` (`response_kind`: `hill`, `linear`, `step`):
//!   - `hill` — uses `half_max` and `hill_power` from columns 5–6; `param_a`/`param_b` ignored.
//!   - `linear` — `param_a` = \(s_\min\), `param_b` = \(s_\max\) for capped linear \(R(s)\).
//!   - `step` — `param_a` = threshold \(s^*\); `param_b` ignored.

use crate::compile::{behavior_rule_set_expr, BehaviorRuleSet};
use crate::response::ResponseFnKind;
use crate::rule::{Response, Rule};

fn parse_dead(s: &str) -> bool {
    matches!(s.trim(), "1" | "true" | "yes")
}

/// Parse one grammar CSV row (8 or 11 columns, semicolon- or comma-separated fields after trim).
pub fn parse_rule_row(line: &str) -> Result<Rule, String> {
    let line = line.trim();
    if line.is_empty() || line.starts_with('#') {
        return Err("empty or comment".into());
    }
    let parts: Vec<&str> = if line.contains(';') {
        line.split(';').map(|p| p.trim()).collect()
    } else {
        line.split(',').map(|p| p.trim()).collect()
    };
    parse_rule_fields(&parts)
}

pub fn parse_rule_fields(parts: &[&str]) -> Result<Rule, String> {
    if parts.len() < 8 {
        return Err(format!("need at least 8 columns, got {}", parts.len()));
    }
    let cell_type = parts[0].to_string();
    let signal = parts[1].to_string();
    let response = match parts[2].to_lowercase().as_str() {
        "increases" => Response::Increases,
        _ => Response::Decreases,
    };
    let behavior = parts[3].to_string();
    let max_response: f64 = parts[4]
        .parse()
        .map_err(|e| format!("max_response: {e}"))?;
    let half_max: f64 = parts[5]
        .parse()
        .map_err(|e| format!("half_max: {e}"))?;
    let hill_power: f64 = parts[6]
        .parse()
        .map_err(|e| format!("hill_power: {e}"))?;
    let applies_to_dead = parse_dead(parts[7]);

    let response_fn = if parts.len() >= 11 {
        match parts[8].to_lowercase().as_str() {
            "linear" => {
                let min_threshold: f64 = parts[9]
                    .parse()
                    .map_err(|e| format!("linear min_threshold: {e}"))?;
                let max_threshold: f64 = parts[10]
                    .parse()
                    .map_err(|e| format!("linear max_threshold: {e}"))?;
                if max_threshold <= min_threshold {
                    return Err("linear requires s_max > s_min".into());
                }
                ResponseFnKind::Linear {
                    min_threshold,
                    max_threshold,
                }
            }
            "step" => {
                let threshold: f64 = parts[9]
                    .parse()
                    .map_err(|e| format!("step threshold: {e}"))?;
                ResponseFnKind::Step { threshold }
            }
            "hill" => {
                if half_max <= 0.0 || hill_power <= 0.0 {
                    return Err("Hill requires positive half_max and hill_power".into());
                }
                ResponseFnKind::Hill {
                    half_max,
                    hill_power,
                }
            }
            other => return Err(format!("unknown response_kind: {other}")),
        }
    } else {
        if half_max <= 0.0 || hill_power <= 0.0 {
            return Err("Hill requires positive half_max and hill_power".into());
        }
        ResponseFnKind::Hill {
            half_max,
            hill_power,
        }
    };

    Ok(Rule {
        cell_type,
        signal,
        response,
        behavior,
        max_response,
        response_fn,
        applies_to_dead,
        condition: None,
    })
}

pub fn rules_from_csv_str(s: &str) -> Result<Vec<Rule>, String> {
    let mut out = Vec::new();
    for line in s.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        match parse_rule_row(line) {
            Ok(r) => out.push(r),
            Err(e) if e == "empty or comment" => {}
            Err(e) => return Err(format!("line {line:?}: {e}")),
        }
    }
    Ok(out)
}

/// Canonical sentence from §3.1 of the *Cell* methods paper (with response details).
pub fn rule_to_grammar_english(r: &Rule) -> String {
    let verb = match r.response {
        Response::Increases => "increases",
        Response::Decreases => "decreases",
    };
    let tail = match &r.response_fn {
        ResponseFnKind::Hill {
            half_max,
            hill_power,
        } => format!(
            " with half-max {} and Hill power {}.",
            half_max, hill_power
        ),
        ResponseFnKind::Linear {
            min_threshold,
            max_threshold,
        } => format!(
            " with linear response from s_min {} to s_max {}.",
            min_threshold, max_threshold
        ),
        ResponseFnKind::Step { threshold } => format!(" with step threshold at {}.", threshold),
        ResponseFnKind::Custom { .. } => ".".to_string(),
    };
    format!(
        "In {}, {} {} {}{}",
        r.cell_type, r.signal, verb, r.behavior, tail
    )
}

/// LaTeX for the bilinear combined behavior \(b(\mathbf{u},\mathbf{d})\) from §5.4 (same [`Expr`](crate::expr::Expr) as evaluation).
pub fn combined_behavior_latex(set: &BehaviorRuleSet) -> String {
    behavior_rule_set_expr(set).simplify().display_latex()
}
