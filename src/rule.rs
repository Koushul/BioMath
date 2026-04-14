use serde::{Deserialize, Serialize};

use crate::expr::Expr;
use crate::response::ResponseFnKind;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Response {
    Increases,
    Decreases,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Rule {
    pub cell_type: String,
    pub signal: String,
    pub response: Response,
    pub behavior: String,
    pub max_response: f64,
    pub response_fn: ResponseFnKind,
    pub applies_to_dead: bool,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub condition: Option<Expr>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct OdeRule {
    pub cell_type: String,
    pub behavior: String,
    pub rate_expr: Expr,
    pub bounds: Option<(f64, f64)>,
}

pub struct RuleBuilder {
    cell_type: String,
    signal: String,
    response: Response,
    behavior: String,
    half_max: Option<f64>,
    hill_power: Option<f64>,
    linear_min: Option<f64>,
    linear_max: Option<f64>,
    step_threshold: Option<f64>,
    max_response: Option<f64>,
    applies_to_dead: bool,
    condition: Option<Expr>,
    custom_response: Option<Expr>,
}

impl Rule {
    pub fn builder(
        cell_type: impl Into<String>,
        signal: impl Into<String>,
        response: Response,
        behavior: impl Into<String>,
    ) -> RuleBuilder {
        RuleBuilder {
            cell_type: cell_type.into(),
            signal: signal.into(),
            response,
            behavior: behavior.into(),
            half_max: None,
            hill_power: None,
            linear_min: None,
            linear_max: None,
            step_threshold: None,
            max_response: None,
            applies_to_dead: false,
            condition: None,
            custom_response: None,
        }
    }

    pub fn to_english(&self) -> String {
        crate::grammar::rule_to_grammar_english(self)
    }

    pub fn to_latex(&self, b_0: f64, b_m: f64, b_min: f64) -> String {
        let s = self.signal.replace('_', "\\_");
        let bname = self.behavior.replace('_', "\\_");
        let r = self.response_fn.to_expr("s").display_latex();
        match self.response {
            Response::Increases => format!(
                "b_{{\\mathrm{{{bname}}}}} = {b0} + ({bM} - {b0}) \\cdot \\left({r}\\right) \\quad (\\mathrm{{with }} s = {s})",
                b0 = b_0,
                bM = b_m,
                r = r.replace("\\mathrm{s}", &format!("\\mathrm{{{s}}}"))
            ),
            Response::Decreases => format!(
                "b_{{\\mathrm{{{bname}}}}} = {b0} - ({b0} - {bm}) \\cdot \\left({r}\\right)",
                b0 = b_0,
                bm = b_min,
                r = r.replace("\\mathrm{s}", &format!("\\mathrm{{{s}}}"))
            ),
        }
    }

    pub fn to_expr(&self, b_0: f64) -> Expr {
        let s = Expr::Signal(self.signal.clone());
        let r = match &self.response_fn {
            ResponseFnKind::Hill {
                half_max,
                hill_power,
            } => Expr::Hill(
                Box::new(s.clone()),
                Box::new(Expr::Const(*half_max)),
                Box::new(Expr::Const(*hill_power)),
            ),
            ResponseFnKind::Linear {
                min_threshold,
                max_threshold,
            } => ResponseFnKind::Linear {
                min_threshold: *min_threshold,
                max_threshold: *max_threshold,
            }
            .to_expr(&self.signal),
            ResponseFnKind::Step { threshold } => ResponseFnKind::Step {
                threshold: *threshold,
            }
            .to_expr(&self.signal),
            ResponseFnKind::Custom { expr } => {
                let mut e = expr.clone();
                e = e.substitute("s", &s);
                e
            }
        };
        let b0 = Expr::Const(b_0);
        match self.response {
            Response::Increases => {
                let bm = Expr::Const(self.max_response);
                b0.clone() + (bm - b0) * r
            }
            Response::Decreases => {
                let bmin = Expr::Const(self.max_response);
                b0.clone() - (b0 - bmin) * r
            }
        }
    }
}

impl OdeRule {
    pub fn to_english(&self) -> String {
        format!(
            "In {}, d {}/dt = {}.",
            self.cell_type,
            self.behavior,
            self.rate_expr.display_math()
        )
    }

    pub fn to_latex(&self) -> String {
        format!(
            "\\frac{{d\\,\\mathrm{{{}}}}}{{dt}} = {}",
            self.behavior.replace('_', "\\_"),
            self.rate_expr.display_latex()
        )
    }
}

impl RuleBuilder {
    pub fn half_max(mut self, v: f64) -> Self {
        self.half_max = Some(v);
        self
    }

    pub fn hill_power(mut self, v: f64) -> Self {
        self.hill_power = Some(v);
        self
    }

    pub fn linear(mut self, min: f64, max: f64) -> Self {
        self.linear_min = Some(min);
        self.linear_max = Some(max);
        self
    }

    pub fn step(mut self, threshold: f64) -> Self {
        self.step_threshold = Some(threshold);
        self
    }

    pub fn max_response(mut self, v: f64) -> Self {
        self.max_response = Some(v);
        self
    }

    pub fn applies_to_dead(mut self, v: bool) -> Self {
        self.applies_to_dead = v;
        self
    }

    pub fn condition(mut self, e: Expr) -> Self {
        self.condition = Some(e);
        self
    }

    pub fn custom_response(mut self, e: Expr) -> Self {
        self.custom_response = Some(e);
        self
    }

    pub fn build(self) -> Rule {
        let response_fn = if let Some(e) = self.custom_response {
            ResponseFnKind::Custom { expr: e }
        } else if let (Some(lo), Some(hi)) = (self.linear_min, self.linear_max) {
            ResponseFnKind::Linear {
                min_threshold: lo,
                max_threshold: hi,
            }
        } else if let Some(t) = self.step_threshold {
            ResponseFnKind::Step { threshold: t }
        } else {
            ResponseFnKind::Hill {
                half_max: self.half_max.expect("half_max required for Hill"),
                hill_power: self.hill_power.expect("hill_power required for Hill"),
            }
        };
        Rule {
            cell_type: self.cell_type,
            signal: self.signal,
            response: self.response,
            behavior: self.behavior,
            max_response: self.max_response.expect("max_response required"),
            response_fn,
            applies_to_dead: self.applies_to_dead,
            condition: self.condition,
        }
    }
}
