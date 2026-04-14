use std::collections::{HashMap, HashSet};

use serde::{Deserialize, Serialize};

use crate::error::CompileError;
use crate::rule::Rule;

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum SignalCategory {
    Chemical,
    Mechanical,
    Contact,
    InternalState,
    Damage,
    Time,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SignalDef {
    pub name: String,
    pub category: SignalCategory,
    pub units: Option<String>,
    pub typical_range: Option<(f64, f64)>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum BehaviorKind {
    Rate,
    Magnitude,
    Probability,
    Transition { target_cell_type: String },
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct BehaviorDef {
    pub name: String,
    pub kind: BehaviorKind,
    pub units: Option<String>,
    pub default_base_value: Option<f64>,
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Dictionary {
    pub signals: HashMap<String, SignalDef>,
    pub behaviors: HashMap<String, BehaviorDef>,
    pub cell_types: HashSet<String>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Warning {
    pub message: String,
}

impl Dictionary {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn register_signal(&mut self, def: SignalDef) {
        self.signals.insert(def.name.clone(), def);
    }

    pub fn register_behavior(&mut self, def: BehaviorDef) {
        self.behaviors.insert(def.name.clone(), def);
    }

    pub fn register_cell_type(&mut self, name: impl Into<String>) {
        self.cell_types.insert(name.into());
    }

    pub fn default_physicell() -> Self {
        let mut d = Dictionary::new();
        for s in [
            ("oxygen", SignalCategory::Chemical),
            ("pressure", SignalCategory::Mechanical),
            ("volume", SignalCategory::InternalState),
            ("ECM", SignalCategory::Mechanical),
            ("cisplatin", SignalCategory::Chemical),
            ("IL-10", SignalCategory::Chemical),
            ("EGF", SignalCategory::Chemical),
            ("time", SignalCategory::Time),
            ("fibroblast_factor", SignalCategory::Chemical),
            ("inflammatory_factor", SignalCategory::Chemical),
            ("pial_contact", SignalCategory::Contact),
        ] {
            d.register_signal(SignalDef {
                name: s.0.into(),
                category: s.1,
                units: None,
                typical_range: None,
            });
        }
        for (name, b0) in [
            ("cycle entry", Some(0.0001_f64)),
            ("necrosis", Some(0.01)),
            ("apoptosis", Some(0.0)),
            ("migration speed", Some(0.0)),
            ("ECM secretion", Some(0.0)),
        ] {
            d.register_behavior(BehaviorDef {
                name: name.into(),
                kind: BehaviorKind::Rate,
                units: None,
                default_base_value: b0,
            });
        }
        for ct in [
            "tumor",
            "motile_tumor",
            "epithelial_tumor",
            "mesenchymal_tumor",
            "fibroblast",
            "MCF-7",
            "naive_T_cell",
            "CD8_T_cell",
            "stem",
            "progenitor",
            "layer_N",
        ] {
            d.register_cell_type(ct);
        }
        d
    }

    pub fn validate_rule(&self, rule: &Rule) -> Result<Vec<Warning>, CompileError> {
        let mut warnings = Vec::new();
        if !self.signals.is_empty() && !self.signals.contains_key(&rule.signal) {
            return Err(CompileError::UnknownSignal {
                name: rule.signal.clone(),
                suggestions: suggestions(&rule.signal, self.signals.keys()),
            });
        }
        if !self.behaviors.is_empty()
            && !self.behaviors.contains_key(&rule.behavior)
            && !is_transition_behavior(&rule.behavior)
        {
            return Err(CompileError::UnknownBehavior {
                name: rule.behavior.clone(),
                suggestions: suggestions(&rule.behavior, self.behaviors.keys()),
            });
        }
        if !self.cell_types.is_empty() && !self.cell_types.contains(&rule.cell_type) {
            return Err(CompileError::UnknownCellType {
                name: rule.cell_type.clone(),
                suggestions: suggestions(&rule.cell_type, self.cell_types.iter()),
            });
        }
        if let Some(BehaviorKind::Transition { target_cell_type }) =
            parse_transition(&rule.behavior)
        {
            if !self.cell_types.is_empty() && !self.cell_types.contains(&target_cell_type) {
                return Err(CompileError::InvalidTransitionTarget {
                    from: rule.cell_type.clone(),
                    to: target_cell_type,
                });
            }
        }
        if self.behaviors.is_empty() && self.signals.is_empty() && self.cell_types.is_empty() {
            return Ok(warnings);
        }
        if let Some(bdef) = self.behaviors.get(&rule.behavior) {
            if bdef.default_base_value.is_none() {
                warnings.push(Warning {
                    message: format!(
                        "no default base value for behavior {}",
                        rule.behavior
                    ),
                });
            }
        }
        Ok(warnings)
    }
}

pub fn is_transition_behavior(behavior: &str) -> bool {
    parse_transition(behavior).is_some()
}

pub fn parse_transition(behavior: &str) -> Option<BehaviorKind> {
    let lower = behavior.to_lowercase();
    let prefixes = ["transition to ", "transformation to "];
    for p in prefixes {
        if let Some(rest) = lower.strip_prefix(p) {
            let target = rest.trim().to_string();
            if !target.is_empty() {
                return Some(BehaviorKind::Transition {
                    target_cell_type: to_snake_case(&target),
                });
            }
        }
    }
    None
}

fn to_snake_case(s: &str) -> String {
    s.chars()
        .map(|c| {
            if c.is_whitespace() || c == '-' {
                '_'
            } else {
                c
            }
        })
        .collect::<String>()
        .split('_')
        .filter(|p| !p.is_empty())
        .collect::<Vec<_>>()
        .join("_")
        .to_lowercase()
}

#[allow(clippy::needless_range_loop)]
fn levenshtein(a: &str, b: &str) -> usize {
    let a: Vec<char> = a.chars().collect();
    let b: Vec<char> = b.chars().collect();
    let n = a.len();
    let m = b.len();
    let mut dp = vec![vec![0usize; m + 1]; n + 1];
    for (i, row) in dp.iter_mut().enumerate().take(n + 1) {
        row[0] = i;
    }
    for j in 1..=m {
        dp[0][j] = j;
    }
    for i in 1..=n {
        for j in 1..=m {
            let cost = usize::from(a[i - 1] != b[j - 1]);
            dp[i][j] = (dp[i - 1][j] + 1)
                .min(dp[i][j - 1] + 1)
                .min(dp[i - 1][j - 1] + cost);
        }
    }
    dp[n][m]
}

fn suggestions<'a>(name: &str, candidates: impl Iterator<Item = &'a String>) -> Vec<String> {
    let mut scored: Vec<(usize, String)> = candidates
        .map(|c| (levenshtein(&name.to_lowercase(), &c.to_lowercase()), c.clone()))
        .filter(|(d, _)| *d <= 5)
        .collect();
    scored.sort_by_key(|(d, _)| *d);
    scored.into_iter().take(5).map(|(_, s)| s).collect()
}

impl Dictionary {
    pub fn suggestions_for_signal(&self, name: &str) -> Vec<String> {
        suggestions(name, self.signals.keys())
    }
}
