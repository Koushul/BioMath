use crate::error::BioMathError;
use crate::expr::EvalContext;
use crate::rule::OdeRule;

pub trait Integrator: Send + Sync {
    fn step(
        &self,
        ode: &OdeRule,
        ctx: &EvalContext,
        current: f64,
        dt: f64,
    ) -> Result<f64, BioMathError>;
}

pub struct Euler;

impl Integrator for Euler {
    fn step(
        &self,
        ode: &OdeRule,
        ctx: &EvalContext,
        current: f64,
        dt: f64,
    ) -> Result<f64, BioMathError> {
        let mut c = ctx.clone();
        c.behaviors.insert(ode.behavior.clone(), current);
        let k = ode.rate_expr.eval(&c)?;
        let mut next = current + dt * k;
        if let Some((lo, hi)) = ode.bounds {
            next = next.clamp(lo, hi);
        }
        Ok(next)
    }
}

pub struct Rk4;

impl Integrator for Rk4 {
    fn step(
        &self,
        ode: &OdeRule,
        ctx: &EvalContext,
        current: f64,
        dt: f64,
    ) -> Result<f64, BioMathError> {
        let f = |y: f64| -> Result<f64, BioMathError> {
            let mut c = ctx.clone();
            c.behaviors.insert(ode.behavior.clone(), y);
            ode.rate_expr.eval(&c)
        };
        let k1 = f(current)?;
        let k2 = f(current + 0.5 * dt * k1)?;
        let k3 = f(current + 0.5 * dt * k2)?;
        let k4 = f(current + dt * k3)?;
        let mut next = current + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        if let Some((lo, hi)) = ode.bounds {
            next = next.clamp(lo, hi);
        }
        Ok(next)
    }
}
