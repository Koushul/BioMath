use std::collections::{HashMap, HashSet};

use serde::{Deserialize, Serialize};

use crate::dictionary::{parse_transition, Dictionary};
use crate::error::CompileError;
use crate::expr::Expr;
use crate::response::{r_from_response, t_from_r, ResponseFnKind};
use crate::rule::{Response, Rule};

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct BaseValueMap(pub HashMap<(String, String), f64>);

impl BaseValueMap {
    pub fn insert(&mut self, cell: impl Into<String>, behavior: impl Into<String>, v: f64) {
        self.0.insert((cell.into(), behavior.into()), v);
    }

    pub fn get(&self, cell: &str, behavior: &str) -> Option<f64> {
        self.0.get(&(cell.to_string(), behavior.to_string())).copied()
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct CompiledUpRule {
    pub signal: String,
    pub response_fn: ResponseFnKind,
    pub max_response: f64,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct CompiledDownRule {
    pub signal: String,
    pub response_fn: ResponseFnKind,
    pub min_response: f64,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct BehaviorRuleSet {
    pub cell_type: String,
    pub behavior: String,
    pub base_value: f64,
    pub max_value: f64,
    pub min_value: f64,
    pub up_rules: Vec<CompiledUpRule>,
    pub down_rules: Vec<CompiledDownRule>,
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct TransitionGraph {
    pub edges: HashMap<String, Vec<String>>,
}

impl TransitionGraph {
    pub fn add(&mut self, from: String, to: String) {
        self.edges.entry(from).or_default().push(to);
    }

    pub fn reachable_from(&self, cell_type: &str) -> HashSet<String> {
        let mut out = HashSet::new();
        let mut stack = vec![cell_type.to_string()];
        while let Some(c) = stack.pop() {
            if !out.insert(c.clone()) {
                continue;
            }
            if let Some(next) = self.edges.get(&c) {
                for t in next {
                    stack.push(t.clone());
                }
            }
        }
        out
    }

    pub fn is_valid_transition(&self, from: &str, to: &str) -> bool {
        self.edges
            .get(from)
            .map(|v| v.iter().any(|x| x == to))
            .unwrap_or(false)
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct CompiledModel {
    pub behavior_sets: Vec<BehaviorRuleSet>,
    pub transitions: TransitionGraph,
    pub raw_rules: Vec<Rule>,
}

pub fn compile(
    rules: Vec<Rule>,
    base_values: &BaseValueMap,
    dictionary: Option<&Dictionary>,
) -> Result<CompiledModel, Vec<CompileError>> {
    let mut errors = Vec::new();
    for (i, rule) in rules.iter().enumerate() {
        if let Some(dict) = dictionary {
            if let Err(e) = dict.validate_rule(rule) {
                errors.push(e);
                continue;
            }
        }
        match &rule.response_fn {
            ResponseFnKind::Hill {
                half_max,
                hill_power,
            } => {
                if *half_max < 0.0 {
                    errors.push(CompileError::NegativeHalfMax { rule_index: i });
                }
                if *half_max == 0.0 {
                    errors.push(CompileError::NonPositiveHalfMax { rule_index: i });
                }
                if *hill_power < 0.0 {
                    errors.push(CompileError::NegativeHillPower { rule_index: i });
                }
            }
            ResponseFnKind::Linear {
                min_threshold,
                max_threshold,
            } => {
                if max_threshold <= min_threshold {
                    errors.push(CompileError::InvalidLinearThresholds { rule_index: i });
                }
            }
            ResponseFnKind::Step { .. } | ResponseFnKind::Custom { .. } => {}
        }
    }
    if !errors.is_empty() {
        return Err(errors);
    }

    let mut groups: HashMap<(String, String), Vec<&Rule>> = HashMap::new();
    for r in &rules {
        groups
            .entry((r.cell_type.clone(), r.behavior.clone()))
            .or_default()
            .push(r);
    }

    let mut behavior_sets = Vec::new();
    let mut transitions = TransitionGraph::default();

    for ((cell_type, behavior), group) in groups {
        let mut up: Vec<&Rule> = Vec::new();
        let mut down: Vec<&Rule> = Vec::new();
        for r in &group {
            match r.response {
                Response::Increases => up.push(r),
                Response::Decreases => down.push(r),
            }
        }

        let mut max_vals: Vec<f64> = up.iter().map(|r| r.max_response).collect();
        max_vals.sort_by(|a, b| a.partial_cmp(b).unwrap());
        max_vals.dedup_by(|a, b| (*a - *b).abs() < 1e-12);
        if max_vals.len() > 1 {
            errors.push(CompileError::InconsistentMaxResponse {
                cell_type: cell_type.clone(),
                behavior: behavior.clone(),
                values: max_vals,
            });
            continue;
        }
        let mut min_vals: Vec<f64> = down.iter().map(|r| r.max_response).collect();
        min_vals.sort_by(|a, b| a.partial_cmp(b).unwrap());
        min_vals.dedup_by(|a, b| (*a - *b).abs() < 1e-12);
        if min_vals.len() > 1 {
            errors.push(CompileError::InconsistentMinResponse {
                cell_type: cell_type.clone(),
                behavior: behavior.clone(),
                values: min_vals,
            });
            continue;
        }

        let b_0 = base_values
            .get(&cell_type, &behavior)
            .or_else(|| {
                dictionary.and_then(|d| {
                    d.behaviors
                        .get(&behavior)
                        .and_then(|b| b.default_base_value)
                })
            })
            .unwrap_or(0.0);

        let b_m = if up.is_empty() {
            b_0
        } else {
            max_vals[0]
        };

        let b_min = if down.is_empty() {
            b_0
        } else {
            min_vals[0]
        };

        if !up.is_empty() && b_0 > b_m + 1e-9 {
            errors.push(CompileError::InvalidBounds {
                cell_type: cell_type.clone(),
                behavior: behavior.clone(),
                base: b_0,
                max: b_m,
                min: b_min,
            });
            continue;
        }
        if !down.is_empty() && b_0 < b_min - 1e-9 {
            errors.push(CompileError::InvalidBounds {
                cell_type: cell_type.clone(),
                behavior: behavior.clone(),
                base: b_0,
                max: b_m,
                min: b_min,
            });
            continue;
        }

        let up_rules: Vec<CompiledUpRule> = up
            .iter()
            .map(|r| CompiledUpRule {
                signal: r.signal.clone(),
                response_fn: r.response_fn.clone(),
                max_response: r.max_response,
            })
            .collect();
        let down_rules: Vec<CompiledDownRule> = down
            .iter()
            .map(|r| CompiledDownRule {
                signal: r.signal.clone(),
                response_fn: r.response_fn.clone(),
                min_response: r.max_response,
            })
            .collect();

        if let Some(crate::dictionary::BehaviorKind::Transition { target_cell_type }) =
            parse_transition(&behavior)
        {
            transitions.add(cell_type.clone(), target_cell_type);
        }

        behavior_sets.push(BehaviorRuleSet {
            cell_type,
            behavior,
            base_value: b_0,
            max_value: b_m,
            min_value: b_min,
            up_rules,
            down_rules,
        });
    }

    if !errors.is_empty() {
        return Err(errors);
    }

    behavior_sets.sort_by(|a, b| {
        a.cell_type
            .cmp(&b.cell_type)
            .then_with(|| a.behavior.cmp(&b.behavior))
    });

    Ok(CompiledModel {
        behavior_sets,
        transitions,
        raw_rules: rules,
    })
}

pub fn pooled_u(signals: &HashMap<String, f64>, rules: &[CompiledUpRule]) -> f64 {
    let mut sum_t = 0.0;
    for r in rules {
        let s = signals.get(&r.signal).copied().unwrap_or(0.0).max(0.0);
        let ri = r_from_response(s, &r.response_fn).unwrap_or(0.0);
        sum_t += t_from_r(ri);
    }
    sum_t / (1.0 + sum_t)
}

pub fn pooled_d(signals: &HashMap<String, f64>, rules: &[CompiledDownRule]) -> f64 {
    let mut sum_t = 0.0;
    for r in rules {
        let s = signals.get(&r.signal).copied().unwrap_or(0.0).max(0.0);
        let ri = r_from_response(s, &r.response_fn).unwrap_or(0.0);
        sum_t += t_from_r(ri);
    }
    sum_t / (1.0 + sum_t)
}

pub fn bilinear(b_0: f64, b_m: f64, b_min: f64, u: f64, d: f64) -> f64 {
    (1.0 - d) * ((1.0 - u) * b_0 + u * b_m) + d * b_min
}

pub fn pooled_u_expr(rules: &[CompiledUpRule]) -> Expr {
    if rules.is_empty() {
        return Expr::Const(0.0);
    }
    let mut sum_t = Expr::Const(0.0);
    for r in rules {
        let ri = r.response_fn.to_expr(&r.signal);
        let denom = Expr::Const(1.0) - ri.clone();
        let t = ri / denom;
        sum_t = sum_t + t;
    }
    sum_t.clone() / (Expr::Const(1.0) + sum_t)
}

pub fn pooled_d_expr(rules: &[CompiledDownRule]) -> Expr {
    if rules.is_empty() {
        return Expr::Const(0.0);
    }
    let mut sum_t = Expr::Const(0.0);
    for r in rules {
        let ri = r.response_fn.to_expr(&r.signal);
        let denom = Expr::Const(1.0) - ri.clone();
        let t = ri / denom;
        sum_t = sum_t + t;
    }
    sum_t.clone() / (Expr::Const(1.0) + sum_t)
}

pub fn behavior_rule_set_expr(set: &BehaviorRuleSet) -> Expr {
    let u = pooled_u_expr(&set.up_rules);
    let d = pooled_d_expr(&set.down_rules);
    let b0 = Expr::Const(set.base_value);
    let bm = Expr::Const(set.max_value);
    let bmin = Expr::Const(set.min_value);
    (Expr::Const(1.0) - d.clone()) * ((Expr::Const(1.0) - u.clone()) * b0 + u * bm) + d * bmin
}
