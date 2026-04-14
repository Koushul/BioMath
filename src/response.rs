use serde::{Deserialize, Serialize};

use crate::error::BioMathError;
use crate::expr::{EvalContext, Expr};

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum ResponseFnKind {
    Hill {
        half_max: f64,
        hill_power: f64,
    },
    Linear {
        min_threshold: f64,
        max_threshold: f64,
    },
    Step {
        threshold: f64,
    },
    Custom {
        expr: Expr,
    },
}

impl ResponseFnKind {
    pub fn evaluate(&self, signal_val: f64) -> Result<f64, BioMathError> {
        r_from_response(signal_val, self)
    }

    pub fn to_expr(&self, signal_var: &str) -> Expr {
        let s = Expr::Signal(signal_var.to_string());
        match self {
            ResponseFnKind::Hill {
                half_max,
                hill_power,
            } => Expr::Hill(
                Box::new(s),
                Box::new(Expr::Const(*half_max)),
                Box::new(Expr::Const(*hill_power)),
            ),
            ResponseFnKind::Linear {
                min_threshold,
                max_threshold,
            } => {
                let smin = Expr::Const(*min_threshold);
                let smax = Expr::Const(*max_threshold);
                let den = smax.clone() - smin.clone();
                Expr::IfGreater(
                    Box::new(s.clone()),
                    Box::new(smax),
                    Box::new(Expr::Const(1.0)),
                    Box::new(Expr::IfGreater(
                        Box::new(smin.clone() - s.clone()),
                        Box::new(Expr::Const(0.0)),
                        Box::new(Expr::Const(0.0)),
                        Box::new((s - smin) / den),
                    )),
                )
            }
            ResponseFnKind::Step { threshold } => Expr::IfGreater(
                Box::new(s),
                Box::new(Expr::Const(*threshold)),
                Box::new(Expr::Const(1.0)),
                Box::new(Expr::Const(0.0)),
            ),
            ResponseFnKind::Custom { expr } => expr.clone(),
        }
    }
}

pub fn r_from_response(signal_val: f64, rf: &ResponseFnKind) -> Result<f64, BioMathError> {
    match rf {
        ResponseFnKind::Hill {
            half_max,
            hill_power,
        } => {
            if *half_max <= 0.0 || *hill_power <= 0.0 {
                return Err(BioMathError::EvalError("invalid Hill params".into()));
            }
            let s = signal_val.max(0.0);
            let sn = s.powf(*hill_power);
            let hn = half_max.powf(*hill_power);
            Ok(sn / (hn + sn))
        }
        ResponseFnKind::Linear {
            min_threshold,
            max_threshold,
        } => {
            if max_threshold <= min_threshold {
                return Err(BioMathError::EvalError("invalid linear thresholds".into()));
            }
            if signal_val < *min_threshold {
                Ok(0.0)
            } else if signal_val > *max_threshold {
                Ok(1.0)
            } else {
                Ok((signal_val - min_threshold) / (max_threshold - min_threshold))
            }
        }
        ResponseFnKind::Step { threshold } => Ok(if signal_val >= *threshold { 1.0 } else { 0.0 }),
        ResponseFnKind::Custom { expr } => {
            let mut ctx = EvalContext::default();
            ctx.signals.insert("s".into(), signal_val);
            expr.eval(&ctx)
        }
    }
}

pub fn t_from_r(r: f64) -> f64 {
    let r = r.clamp(0.0, 1.0 - 1e-15);
    r / (1.0 - r).max(1e-15)
}
