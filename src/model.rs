use std::collections::HashMap;

use serde::{Deserialize, Serialize};

use crate::compile::{
    behavior_rule_set_expr, bilinear, compile, BaseValueMap, BehaviorRuleSet, CompiledModel,
    pooled_d, pooled_u,
};
use crate::dictionary::Dictionary;
use crate::error::BioMathError;
use crate::expr::{EvalContext, Expr};
use crate::integrate::{Euler, Integrator, Rk4};
use crate::rule::{OdeRule, Rule};

pub type SignalEnv = HashMap<String, f64>;
pub type BehaviorMap = HashMap<String, f64>;

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct CellState {
    pub signals: SignalEnv,
    pub behaviors: BehaviorMap,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Model {
    pub name: String,
    pub compiled: CompiledModel,
    #[serde(default)]
    pub ode_rules: Vec<OdeRule>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct EvalRequest {
    pub cell_type: String,
    pub signals: SignalEnv,
    #[serde(default)]
    pub state: BehaviorMap,
    #[serde(default)]
    pub dt: Option<f64>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct EvalResponse {
    pub behaviors: BehaviorMap,
    #[serde(default, skip_serializing_if = "HashMap::is_empty")]
    pub updated_state: BehaviorMap,
}

impl Model {
    pub fn from_rules(
        name: impl Into<String>,
        rules: Vec<Rule>,
        base_values: BaseValueMap,
        dictionary: Option<&Dictionary>,
    ) -> Result<Self, Vec<crate::error::CompileError>> {
        let compiled = compile(rules, &base_values, dictionary)?;
        Ok(Model {
            name: name.into(),
            compiled,
            ode_rules: Vec::new(),
        })
    }

    pub fn with_ode_rules(mut self, odes: Vec<OdeRule>) -> Self {
        self.ode_rules = odes;
        self
    }

    pub fn evaluate(&self, cell_type: &str, signals: &SignalEnv) -> BehaviorMap {
        let mut out = BehaviorMap::new();
        for set in &self.compiled.behavior_sets {
            if set.cell_type != cell_type {
                continue;
            }
            let u = pooled_u(signals, &set.up_rules);
            let d = pooled_d(signals, &set.down_rules);
            let v = bilinear(set.base_value, set.max_value, set.min_value, u, d);
            out.insert(set.behavior.clone(), v);
        }
        out
    }

    pub fn evaluate_with_conditions(
        &self,
        cell_type: &str,
        signals: &SignalEnv,
    ) -> BehaviorMap {
        let base = self.evaluate(cell_type, signals);
        for r in &self.compiled.raw_rules {
            if r.cell_type != cell_type {
                continue;
            }
            if let Some(cond) = &r.condition {
                let mut ctx = EvalContext::default();
                for (k, v) in signals {
                    ctx.signals.insert(k.clone(), *v);
                }
                if cond.eval(&ctx).unwrap_or(0.0) <= 0.0 {
                    continue;
                }
            }
        }
        base
    }

    pub fn step(
        &self,
        cell_type: &str,
        state: &mut CellState,
        dt: f64,
        integrator: &dyn Integrator,
    ) -> Result<(), BioMathError> {
        let mut ctx = EvalContext::default();
        for (k, v) in &state.signals {
            ctx.signals.insert(k.clone(), *v);
        }
        for (k, v) in &state.behaviors {
            ctx.behaviors.insert(k.clone(), *v);
        }
        for ode in &self.ode_rules {
            if ode.cell_type != cell_type {
                continue;
            }
            let cur = state.behaviors.get(&ode.behavior).copied().unwrap_or(0.0);
            let next = integrator.step(ode, &ctx, cur, dt)?;
            state.behaviors.insert(ode.behavior.clone(), next);
            ctx.behaviors.insert(ode.behavior.clone(), next);
        }
        Ok(())
    }

    pub fn step_euler(&self, cell_type: &str, state: &mut CellState, dt: f64) -> Result<(), BioMathError> {
        self.step(cell_type, state, dt, &Euler)
    }

    pub fn step_rk4(&self, cell_type: &str, state: &mut CellState, dt: f64) -> Result<(), BioMathError> {
        self.step(cell_type, state, dt, &Rk4)
    }

    pub fn behavior_rule_set(&self, cell_type: &str, behavior: &str) -> Option<&BehaviorRuleSet> {
        self.compiled
            .behavior_sets
            .iter()
            .find(|s| s.cell_type == cell_type && s.behavior == behavior)
    }

    pub fn sensitivity(&self, cell_type: &str, behavior: &str, wrt: &str) -> Option<Expr> {
        let set = self.behavior_rule_set(cell_type, behavior)?;
        Some(behavior_rule_set_expr(set).diff(wrt).simplify())
    }

    pub fn sensitivity_at(
        &self,
        cell_type: &str,
        behavior: &str,
        wrt: &str,
        env: &SignalEnv,
        eps: f64,
    ) -> Option<f64> {
        let mut e1 = env.clone();
        let mut e2 = env.clone();
        let v = *e1.get(wrt)?;
        e1.insert(wrt.to_string(), v + eps);
        e2.insert(wrt.to_string(), v - eps);
        let b1 = *self.evaluate(cell_type, &e1).get(behavior)?;
        let b2 = *self.evaluate(cell_type, &e2).get(behavior)?;
        Some((b1 - b2) / (2.0 * eps))
    }

    pub fn jacobian(&self, cell_type: &str) -> Vec<Vec<Expr>> {
        let mut behs: Vec<String> = self
            .ode_rules
            .iter()
            .filter(|o| o.cell_type == cell_type)
            .map(|o| o.behavior.clone())
            .collect();
        behs.sort();
        behs.dedup();
        let n = behs.len();
        let mut m = vec![vec![Expr::Const(0.0); n]; n];
        for (i, bi) in behs.iter().enumerate() {
            if let Some(ode) = self
                .ode_rules
                .iter()
                .find(|o| o.cell_type == cell_type && o.behavior == *bi)
            {
                for (j, bj) in behs.iter().enumerate() {
                    m[i][j] = ode.rate_expr.diff(bj).simplify();
                }
            }
        }
        m
    }

    pub fn jacobian_at(
        &self,
        cell_type: &str,
        state: &CellState,
    ) -> Result<Vec<Vec<f64>>, BioMathError> {
        let j = self.jacobian(cell_type);
        let mut ctx = EvalContext::default();
        for (k, v) in &state.signals {
            ctx.signals.insert(k.clone(), *v);
        }
        for (k, v) in &state.behaviors {
            ctx.behaviors.insert(k.clone(), *v);
        }
        j.into_iter()
            .map(|row| {
                row.into_iter()
                    .map(|e| e.eval(&ctx))
                    .collect::<Result<Vec<_>, _>>()
            })
            .collect::<Result<Vec<_>, _>>()
    }

    pub fn to_json(&self) -> Result<String, BioMathError> {
        serde_json::to_string(self).map_err(|e| BioMathError::Json(e.to_string()))
    }

    pub fn from_json(s: &str) -> Result<Self, BioMathError> {
        serde_json::from_str(s).map_err(|e| BioMathError::Json(e.to_string()))
    }
}
