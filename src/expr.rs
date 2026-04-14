use std::collections::HashSet;
use std::ops::{Add, Div, Mul, Neg, Sub};

use serde::{Deserialize, Serialize};

use crate::error::BioMathError;

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum Expr {
    Const(f64),
    Signal(String),
    Behavior(String),
    Param(String),
    Add(Box<Expr>, Box<Expr>),
    Sub(Box<Expr>, Box<Expr>),
    Mul(Box<Expr>, Box<Expr>),
    Div(Box<Expr>, Box<Expr>),
    Neg(Box<Expr>),
    Pow(Box<Expr>, Box<Expr>),
    Hill(Box<Expr>, Box<Expr>, Box<Expr>),
    Exp(Box<Expr>),
    Log(Box<Expr>),
    Abs(Box<Expr>),
    Min(Box<Expr>, Box<Expr>),
    Max(Box<Expr>, Box<Expr>),
    Clamp(Box<Expr>, Box<Expr>, Box<Expr>),
    IfGreater(Box<Expr>, Box<Expr>, Box<Expr>, Box<Expr>),
}

#[derive(Clone, Debug, Default)]
pub struct EvalContext {
    pub signals: std::collections::HashMap<String, f64>,
    pub behaviors: std::collections::HashMap<String, f64>,
    pub params: std::collections::HashMap<String, f64>,
}

impl Expr {
    pub fn eval(&self, ctx: &EvalContext) -> Result<f64, BioMathError> {
        match self {
            Expr::Const(c) => Ok(*c),
            Expr::Signal(name) => ctx
                .signals
                .get(name)
                .copied()
                .ok_or_else(|| BioMathError::UnknownSignal {
                    name: name.clone(),
                }),
            Expr::Behavior(name) => ctx
                .behaviors
                .get(name)
                .copied()
                .ok_or_else(|| BioMathError::UnknownBehavior {
                    name: name.clone(),
                }),
            Expr::Param(name) => ctx
                .params
                .get(name)
                .copied()
                .ok_or_else(|| BioMathError::UnknownParam {
                    name: name.clone(),
                }),
            Expr::Add(a, b) => Ok(a.eval(ctx)? + b.eval(ctx)?),
            Expr::Sub(a, b) => Ok(a.eval(ctx)? - b.eval(ctx)?),
            Expr::Mul(a, b) => Ok(a.eval(ctx)? * b.eval(ctx)?),
            Expr::Div(a, b) => {
                let den = b.eval(ctx)?;
                if den == 0.0 {
                    return Err(BioMathError::DivisionByZero);
                }
                Ok(a.eval(ctx)? / den)
            }
            Expr::Neg(a) => Ok(-a.eval(ctx)?),
            Expr::Pow(a, b) => {
                let x = a.eval(ctx)?;
                let y = b.eval(ctx)?;
                Ok(x.powf(y))
            }
            Expr::Hill(s, h, n) => {
                let sv = s.eval(ctx)?;
                let hv = h.eval(ctx)?;
                let nv = n.eval(ctx)?;
                if hv <= 0.0 || nv <= 0.0 {
                    return Err(BioMathError::EvalError(
                        "Hill requires positive half_max and hill_power".into(),
                    ));
                }
                let sn = sv.powf(nv);
                let hn = hv.powf(nv);
                Ok(sn / (hn + sn))
            }
            Expr::Exp(a) => Ok(a.eval(ctx)?.exp()),
            Expr::Log(a) => {
                let x = a.eval(ctx)?;
                if x <= 0.0 {
                    return Err(BioMathError::LogDomain);
                }
                Ok(x.ln())
            }
            Expr::Abs(a) => Ok(a.eval(ctx)?.abs()),
            Expr::Min(a, b) => Ok(a.eval(ctx)?.min(b.eval(ctx)?)),
            Expr::Max(a, b) => Ok(a.eval(ctx)?.max(b.eval(ctx)?)),
            Expr::Clamp(v, lo, hi) => {
                let x = v.eval(ctx)?;
                let l = lo.eval(ctx)?;
                let h = hi.eval(ctx)?;
                Ok(x.clamp(l, h))
            }
            Expr::IfGreater(a, b, c, d) => {
                if a.eval(ctx)? > b.eval(ctx)? {
                    c.eval(ctx)
                } else {
                    d.eval(ctx)
                }
            }
        }
    }

    pub fn diff(&self, var: &str) -> Expr {
        match self {
            Expr::Const(_) => Expr::Const(0.0),
            Expr::Signal(n) | Expr::Behavior(n) | Expr::Param(n) => {
                if n == var {
                    Expr::Const(1.0)
                } else {
                    Expr::Const(0.0)
                }
            }
            Expr::Add(a, b) => a.diff(var) + b.diff(var),
            Expr::Sub(a, b) => a.diff(var) - b.diff(var),
            Expr::Mul(a, b) => {
                a.diff(var).mul_ref(b.as_ref()) + a.as_ref().mul_ref(&b.diff(var))
            }
            Expr::Div(a, b) => {
                let num = a.diff(var).mul_ref(b.as_ref()) - a.as_ref().mul_ref(&b.diff(var));
                num / b.as_ref().mul_ref(b.as_ref())
            }
            Expr::Neg(a) => -a.as_ref().diff(var),
            Expr::Pow(base, exp) => Expr::diff_pow(base.as_ref(), exp.as_ref(), var),
            Expr::Hill(s, h, n) => Expr::hill_primitive(s.as_ref(), h.as_ref(), n.as_ref()).diff(var),
            Expr::Exp(x) => Expr::Exp(x.clone()) * x.as_ref().diff(var),
            Expr::Log(x) => x.as_ref().diff(var) / (**x).clone(),
            Expr::Abs(x) => Expr::IfGreater(
                Box::new((**x).clone()),
                Box::new(Expr::Const(0.0)),
                Box::new(x.as_ref().diff(var)),
                Box::new(-x.as_ref().diff(var)),
            ),
            Expr::Min(a, b) => Expr::IfGreater(
                Box::new((**a).clone() - (**b).clone()),
                Box::new(Expr::Const(0.0)),
                Box::new(b.as_ref().diff(var)),
                Box::new(a.as_ref().diff(var)),
            ),
            Expr::Max(a, b) => Expr::IfGreater(
                Box::new((**a).clone() - (**b).clone()),
                Box::new(Expr::Const(0.0)),
                Box::new(a.as_ref().diff(var)),
                Box::new(b.as_ref().diff(var)),
            ),
            Expr::Clamp(v, lo, hi) => {
                let dv = v.as_ref().diff(var);
                let in_range = Expr::IfGreater(
                    Box::new((**v).clone() - (**lo).clone()),
                    Box::new(Expr::Const(0.0)),
                    Box::new(Expr::Const(1.0)),
                    Box::new(Expr::Const(0.0)),
                ) * Expr::IfGreater(
                    Box::new((**hi).clone() - (**v).clone()),
                    Box::new(Expr::Const(0.0)),
                    Box::new(Expr::Const(1.0)),
                    Box::new(Expr::Const(0.0)),
                );
                in_range * dv
            }
            Expr::IfGreater(a, b, c, d) => Expr::IfGreater(
                a.clone(),
                b.clone(),
                Box::new(c.diff(var)),
                Box::new(d.diff(var)),
            ),
        }
    }

    fn diff_pow(base: &Expr, exp: &Expr, var: &str) -> Expr {
        if let Expr::Const(c) = exp {
            return Expr::Const(*c) * base.powf(c - 1.0) * base.diff(var);
        }
        if let Expr::Const(c) = base {
            if *c <= 0.0 {
                return Expr::Const(0.0);
            }
            return Expr::Pow(Box::new(Expr::Const(*c)), Box::new(exp.clone()))
                * Expr::ln_const(*c)
                * exp.diff(var);
        }
        let a = base.clone();
        let b = exp.clone();
        Expr::Pow(Box::new(a.clone()), Box::new(b.clone()))
            * (b.diff(var) * Expr::Log(Box::new(a.clone()))
                + b.clone() * (a.diff(var) / a.clone()))
    }

    fn hill_primitive(s: &Expr, h: &Expr, n: &Expr) -> Expr {
        let num = Expr::Pow(Box::new(s.clone()), Box::new(n.clone()));
        let den = Expr::Pow(Box::new(h.clone()), Box::new(n.clone())) + num.clone();
        num / den
    }

    pub fn expand_hill(&self) -> Expr {
        match self {
            Expr::Hill(s, h, n) => Expr::hill_primitive(s, h, n),
            Expr::Add(a, b) => Expr::Add(
                Box::new(a.expand_hill()),
                Box::new(b.expand_hill()),
            ),
            Expr::Sub(a, b) => Expr::Sub(
                Box::new(a.expand_hill()),
                Box::new(b.expand_hill()),
            ),
            Expr::Mul(a, b) => Expr::Mul(
                Box::new(a.expand_hill()),
                Box::new(b.expand_hill()),
            ),
            Expr::Div(a, b) => Expr::Div(
                Box::new(a.expand_hill()),
                Box::new(b.expand_hill()),
            ),
            Expr::Neg(a) => Expr::Neg(Box::new(a.expand_hill())),
            Expr::Pow(a, b) => Expr::Pow(
                Box::new(a.expand_hill()),
                Box::new(b.expand_hill()),
            ),
            Expr::Exp(a) => Expr::Exp(Box::new(a.expand_hill())),
            Expr::Log(a) => Expr::Log(Box::new(a.expand_hill())),
            Expr::Abs(a) => Expr::Abs(Box::new(a.expand_hill())),
            Expr::Min(a, b) => Expr::Min(
                Box::new(a.expand_hill()),
                Box::new(b.expand_hill()),
            ),
            Expr::Max(a, b) => Expr::Max(
                Box::new(a.expand_hill()),
                Box::new(b.expand_hill()),
            ),
            Expr::Clamp(a, b, c) => Expr::Clamp(
                Box::new(a.expand_hill()),
                Box::new(b.expand_hill()),
                Box::new(c.expand_hill()),
            ),
            Expr::IfGreater(a, b, c, d) => Expr::IfGreater(
                Box::new(a.expand_hill()),
                Box::new(b.expand_hill()),
                Box::new(c.expand_hill()),
                Box::new(d.expand_hill()),
            ),
            _ => self.clone(),
        }
    }

    pub fn simplify(&self) -> Expr {
        let e = match self {
            Expr::Add(a, b) => {
                let a = a.simplify();
                let b = b.simplify();
                match (&a, &b) {
                    (Expr::Const(x), Expr::Const(y)) => Expr::Const(x + y),
                    (Expr::Const(c), _) if *c == 0.0 => b,
                    (_, Expr::Const(c)) if *c == 0.0 => a,
                    _ => Expr::Add(Box::new(a), Box::new(b)),
                }
            }
            Expr::Sub(a, b) => {
                let a = a.simplify();
                let b = b.simplify();
                match (&a, &b) {
                    (Expr::Const(x), Expr::Const(y)) => Expr::Const(x - y),
                    (_, Expr::Const(c)) if *c == 0.0 => a,
                    _ if a.structural_eq(&b) => Expr::Const(0.0),
                    _ => Expr::Sub(Box::new(a), Box::new(b)),
                }
            }
            Expr::Mul(a, b) => {
                let a = a.simplify();
                let b = b.simplify();
                match (&a, &b) {
                    (Expr::Const(x), Expr::Const(y)) => Expr::Const(x * y),
                    (Expr::Const(c), _) | (_, Expr::Const(c)) if *c == 0.0 => Expr::Const(0.0),
                    (Expr::Const(c), _x) | (_x, Expr::Const(c)) if *c == 1.0 => {
                        if matches!(a, Expr::Const(_)) {
                            b
                        } else {
                            a
                        }
                    }
                    _ => Expr::Mul(Box::new(a), Box::new(b)),
                }
            }
            Expr::Div(a, b) => {
                let a = a.simplify();
                let b = b.simplify();
                match (&a, &b) {
                    (Expr::Const(x), Expr::Const(y)) if *y != 0.0 => Expr::Const(x / y),
                    (_, Expr::Const(c)) if *c == 0.0 => Expr::Div(Box::new(a), Box::new(b)),
                    _ if a.structural_eq(&b) => Expr::Const(1.0),
                    (_, Expr::Const(c)) if *c == 1.0 => a,
                    _ => Expr::Div(Box::new(a), Box::new(b)),
                }
            }
            Expr::Neg(a) => {
                let a = a.simplify();
                match a {
                    Expr::Const(c) => Expr::Const(-c),
                    Expr::Neg(x) => *x,
                    _ => Expr::Neg(Box::new(a)),
                }
            }
            Expr::Pow(a, b) => {
                let a = a.simplify();
                let b = b.simplify();
                match (&a, &b) {
                    (Expr::Const(x), Expr::Const(y)) => Expr::Const(x.powf(*y)),
                    (_, Expr::Const(e)) if *e == 0.0 => Expr::Const(1.0),
                    (x, Expr::Const(e)) if *e == 1.0 => x.clone(),
                    _ => Expr::Pow(Box::new(a), Box::new(b)),
                }
            }
            Expr::Exp(a) => {
                let a = a.simplify();
                match a {
                    Expr::Log(x) => *x,
                    Expr::Const(c) => Expr::Const(c.exp()),
                    _ => Expr::Exp(Box::new(a)),
                }
            }
            Expr::Log(a) => {
                let a = a.simplify();
                match a {
                    Expr::Exp(x) => *x,
                    Expr::Const(c) => Expr::Const(c.ln()),
                    _ => Expr::Log(Box::new(a)),
                }
            }
            Expr::Hill(s, h, n) => Expr::Hill(
                Box::new(s.simplify()),
                Box::new(h.simplify()),
                Box::new(n.simplify()),
            ),
            Expr::Abs(a) => Expr::Abs(Box::new(a.simplify())),
            Expr::Min(a, b) => Expr::Min(
                Box::new(a.simplify()),
                Box::new(b.simplify()),
            ),
            Expr::Max(a, b) => Expr::Max(
                Box::new(a.simplify()),
                Box::new(b.simplify()),
            ),
            Expr::Clamp(a, b, c) => Expr::Clamp(
                Box::new(a.simplify()),
                Box::new(b.simplify()),
                Box::new(c.simplify()),
            ),
            Expr::IfGreater(a, b, c, d) => Expr::IfGreater(
                Box::new(a.simplify()),
                Box::new(b.simplify()),
                Box::new(c.simplify()),
                Box::new(d.simplify()),
            ),
            _ => self.clone(),
        };
        if !e.structural_eq(self) {
            e.simplify()
        } else {
            e
        }
    }

    pub fn diff_n(&self, var: &str, n: usize) -> Expr {
        let mut cur = self.clone();
        for _ in 0..n {
            cur = cur.diff(var).simplify();
        }
        cur
    }

    pub fn gradient(&self, vars: &[&str]) -> Vec<Expr> {
        vars.iter().map(|v| self.diff(v).simplify()).collect()
    }

    pub fn free_variables(&self) -> HashSet<String> {
        let mut s = HashSet::new();
        self.collect_vars(&mut s);
        s
    }

    fn collect_vars(&self, out: &mut HashSet<String>) {
        match self {
            Expr::Signal(n) | Expr::Behavior(n) | Expr::Param(n) => {
                out.insert(n.clone());
            }
            Expr::Add(a, b)
            | Expr::Sub(a, b)
            | Expr::Mul(a, b)
            | Expr::Div(a, b)
            | Expr::Pow(a, b)
            | Expr::Min(a, b)
            | Expr::Max(a, b) => {
                a.collect_vars(out);
                b.collect_vars(out);
            }
            Expr::Neg(a)
            | Expr::Exp(a)
            | Expr::Log(a)
            | Expr::Abs(a) => a.collect_vars(out),
            Expr::Hill(x, y, z) => {
                x.collect_vars(out);
                y.collect_vars(out);
                z.collect_vars(out);
            }
            Expr::Clamp(a, b, c) => {
                a.collect_vars(out);
                b.collect_vars(out);
                c.collect_vars(out);
            }
            Expr::IfGreater(a, b, c, d) => {
                a.collect_vars(out);
                b.collect_vars(out);
                c.collect_vars(out);
                d.collect_vars(out);
            }
            Expr::Const(_) => {}
        }
    }

    pub fn substitute(&self, var: &str, replacement: &Expr) -> Expr {
        match self {
            Expr::Signal(n) | Expr::Behavior(n) | Expr::Param(n) if n == var => {
                replacement.clone()
            }
            Expr::Add(a, b) => {
                a.substitute(var, replacement) + b.substitute(var, replacement)
            }
            Expr::Sub(a, b) => {
                a.substitute(var, replacement) - b.substitute(var, replacement)
            }
            Expr::Mul(a, b) => {
                a.substitute(var, replacement).mul_ref(&b.substitute(var, replacement))
            }
            Expr::Div(a, b) => {
                a.substitute(var, replacement) / b.substitute(var, replacement)
            }
            Expr::Neg(a) => -a.substitute(var, replacement),
            Expr::Pow(a, b) => Expr::Pow(
                Box::new(a.substitute(var, replacement)),
                Box::new(b.substitute(var, replacement)),
            ),
            Expr::Hill(a, b, c) => Expr::Hill(
                Box::new(a.substitute(var, replacement)),
                Box::new(b.substitute(var, replacement)),
                Box::new(c.substitute(var, replacement)),
            ),
            Expr::Exp(a) => Expr::Exp(Box::new(a.substitute(var, replacement))),
            Expr::Log(a) => Expr::Log(Box::new(a.substitute(var, replacement))),
            Expr::Abs(a) => Expr::Abs(Box::new(a.substitute(var, replacement))),
            Expr::Min(a, b) => Expr::Min(
                Box::new(a.substitute(var, replacement)),
                Box::new(b.substitute(var, replacement)),
            ),
            Expr::Max(a, b) => Expr::Max(
                Box::new(a.substitute(var, replacement)),
                Box::new(b.substitute(var, replacement)),
            ),
            Expr::Clamp(a, b, c) => Expr::Clamp(
                Box::new(a.substitute(var, replacement)),
                Box::new(b.substitute(var, replacement)),
                Box::new(c.substitute(var, replacement)),
            ),
            Expr::IfGreater(a, b, c, d) => Expr::IfGreater(
                Box::new(a.substitute(var, replacement)),
                Box::new(b.substitute(var, replacement)),
                Box::new(c.substitute(var, replacement)),
                Box::new(d.substitute(var, replacement)),
            ),
            _ => self.clone(),
        }
    }

    pub fn is_const(&self) -> bool {
        matches!(self, Expr::Const(_))
    }

    pub fn structural_eq(&self, other: &Expr) -> bool {
        match (self, other) {
            (Expr::Const(a), Expr::Const(b)) => a == b,
            (Expr::Signal(x), Expr::Signal(y)) => x == y,
            (Expr::Behavior(x), Expr::Behavior(y)) => x == y,
            (Expr::Param(x), Expr::Param(y)) => x == y,
            (Expr::Add(a1, b1), Expr::Add(a2, b2)) => {
                a1.structural_eq(a2) && b1.structural_eq(b2)
            }
            (Expr::Sub(a1, b1), Expr::Sub(a2, b2)) => {
                a1.structural_eq(a2) && b1.structural_eq(b2)
            }
            (Expr::Mul(a1, b1), Expr::Mul(a2, b2)) => {
                a1.structural_eq(a2) && b1.structural_eq(b2)
            }
            (Expr::Div(a1, b1), Expr::Div(a2, b2)) => {
                a1.structural_eq(a2) && b1.structural_eq(b2)
            }
            (Expr::Neg(a1), Expr::Neg(a2)) => a1.structural_eq(a2),
            (Expr::Pow(a1, b1), Expr::Pow(a2, b2)) => {
                a1.structural_eq(a2) && b1.structural_eq(b2)
            }
            (Expr::Hill(a1, b1, c1), Expr::Hill(a2, b2, c2)) => {
                a1.structural_eq(a2) && b1.structural_eq(b2) && c1.structural_eq(c2)
            }
            (Expr::Exp(a1), Expr::Exp(a2)) => a1.structural_eq(a2),
            (Expr::Log(a1), Expr::Log(a2)) => a1.structural_eq(a2),
            (Expr::Abs(a1), Expr::Abs(a2)) => a1.structural_eq(a2),
            (Expr::Min(a1, b1), Expr::Min(a2, b2)) => {
                a1.structural_eq(a2) && b1.structural_eq(b2)
            }
            (Expr::Max(a1, b1), Expr::Max(a2, b2)) => {
                a1.structural_eq(a2) && b1.structural_eq(b2)
            }
            (Expr::Clamp(a1, b1, c1), Expr::Clamp(a2, b2, c2)) => {
                a1.structural_eq(a2) && b1.structural_eq(b2) && c1.structural_eq(c2)
            }
            (Expr::IfGreater(a1, b1, c1, d1), Expr::IfGreater(a2, b2, c2, d2)) => {
                a1.structural_eq(a2)
                    && b1.structural_eq(b2)
                    && c1.structural_eq(c2)
                    && d1.structural_eq(d2)
            }
            _ => false,
        }
    }

    pub fn display_math(&self) -> String {
        match self {
            Expr::Const(c) => format!("{c}"),
            Expr::Signal(n) | Expr::Behavior(n) | Expr::Param(n) => n.clone(),
            Expr::Add(a, b) => format!("({} + {})", a.display_math(), b.display_math()),
            Expr::Sub(a, b) => format!("({} - {})", a.display_math(), b.display_math()),
            Expr::Mul(a, b) => format!("({} * {})", a.display_math(), b.display_math()),
            Expr::Div(a, b) => format!("({} / {})", a.display_math(), b.display_math()),
            Expr::Neg(a) => format!("(-{})", a.display_math()),
            Expr::Pow(a, b) => format!("({} ^ {})", a.display_math(), b.display_math()),
            Expr::Hill(s, h, n) => format!(
                "Hill({}, {}, {})",
                s.display_math(),
                h.display_math(),
                n.display_math()
            ),
            Expr::Exp(a) => format!("exp({})", a.display_math()),
            Expr::Log(a) => format!("log({})", a.display_math()),
            Expr::Abs(a) => format!("abs({})", a.display_math()),
            Expr::Min(a, b) => format!("min({}, {})", a.display_math(), b.display_math()),
            Expr::Max(a, b) => format!("max({}, {})", a.display_math(), b.display_math()),
            Expr::Clamp(v, lo, hi) => format!(
                "clamp({}, {}, {})",
                v.display_math(),
                lo.display_math(),
                hi.display_math()
            ),
            Expr::IfGreater(a, b, c, d) => format!(
                "if {} > {} then {} else {}",
                a.display_math(),
                b.display_math(),
                c.display_math(),
                d.display_math()
            ),
        }
    }

    pub fn display_latex(&self) -> String {
        match self {
            Expr::Const(c) => format!("{c}"),
            Expr::Signal(n) | Expr::Behavior(n) | Expr::Param(n) => {
                format!("\\mathrm{{{}}}", n.replace('_', "\\_"))
            }
            Expr::Add(a, b) => format!("{} + {}", a.display_latex(), b.display_latex()),
            Expr::Sub(a, b) => format!("{} - {}", a.display_latex(), b.display_latex()),
            Expr::Mul(a, b) => format!("{} \\cdot {}", a.display_latex(), b.display_latex()),
            Expr::Div(a, b) => format!(
                "\\frac{{{}}}{{{}}}",
                a.display_latex(),
                b.display_latex()
            ),
            Expr::Neg(a) => format!("-{}", a.display_latex()),
            Expr::Pow(a, b) => format!("{}^{{{}}}", a.display_latex(), b.display_latex()),
            Expr::Hill(s, h, n) => format!(
                "\\frac{{{}^{{{}}}}}{{{}^{{{}}} + {}^{{{}}}}}",
                s.display_latex(),
                n.display_latex(),
                h.display_latex(),
                n.display_latex(),
                s.display_latex(),
                n.display_latex(),
            ),
            Expr::Exp(a) => format!("\\exp({})", a.display_latex()),
            Expr::Log(a) => format!("\\ln({})", a.display_latex()),
            Expr::Abs(a) => format!("\\left|{}\\right|", a.display_latex()),
            Expr::Min(a, b) => format!(
                "\\min({}, {})",
                a.display_latex(),
                b.display_latex()
            ),
            Expr::Max(a, b) => format!(
                "\\max({}, {})",
                a.display_latex(),
                b.display_latex()
            ),
            Expr::Clamp(v, lo, hi) => format!(
                "\\mathrm{{clamp}}({}, {}, {})",
                v.display_latex(),
                lo.display_latex(),
                hi.display_latex()
            ),
            Expr::IfGreater(a, b, c, d) => format!(
                "\\begin{{cases}} {} & {} > {} \\\\ {} & \\mathrm{{else}} \\end{{cases}}",
                c.display_latex(),
                a.display_latex(),
                b.display_latex(),
                d.display_latex()
            ),
        }
    }

    fn mul_ref(&self, other: &Expr) -> Expr {
        Expr::Mul(Box::new(self.clone()), Box::new(other.clone()))
    }

    fn powf(&self, exp: f64) -> Expr {
        Expr::Pow(Box::new(self.clone()), Box::new(Expr::Const(exp)))
    }

    fn ln_const(c: f64) -> Expr {
        Expr::Const(c.ln())
    }
}

pub fn const_(c: f64) -> Expr {
    Expr::Const(c)
}

pub fn signal(name: impl Into<String>) -> Expr {
    Expr::Signal(name.into())
}

pub fn behavior(name: impl Into<String>) -> Expr {
    Expr::Behavior(name.into())
}

pub fn param(name: impl Into<String>) -> Expr {
    Expr::Param(name.into())
}

pub fn hill(signal_expr: Expr, half_max: f64, power: f64) -> Expr {
    Expr::Hill(
        Box::new(signal_expr),
        Box::new(Expr::Const(half_max)),
        Box::new(Expr::Const(power)),
    )
}

impl Add for Expr {
    type Output = Expr;
    fn add(self, rhs: Expr) -> Expr {
        Expr::Add(Box::new(self), Box::new(rhs))
    }
}

impl Sub for Expr {
    type Output = Expr;
    fn sub(self, rhs: Expr) -> Expr {
        Expr::Sub(Box::new(self), Box::new(rhs))
    }
}

impl Mul for Expr {
    type Output = Expr;
    fn mul(self, rhs: Expr) -> Expr {
        Expr::Mul(Box::new(self), Box::new(rhs))
    }
}

impl Div for Expr {
    type Output = Expr;
    fn div(self, rhs: Expr) -> Expr {
        Expr::Div(Box::new(self), Box::new(rhs))
    }
}

impl Neg for Expr {
    type Output = Expr;
    fn neg(self) -> Expr {
        Expr::Neg(Box::new(self))
    }
}
