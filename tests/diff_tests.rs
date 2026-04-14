use approx::assert_relative_eq;
use bio_math::expr::{const_, hill, signal, EvalContext, Expr};

fn ctx_s(name: &str, v: f64) -> EvalContext {
    let mut c = EvalContext::default();
    c.signals.insert(name.into(), v);
    c
}

fn ctx_multi(pairs: &[(&str, f64)]) -> EvalContext {
    let mut c = EvalContext::default();
    for (k, v) in pairs {
        c.signals.insert((*k).into(), *v);
    }
    c
}

fn verify_fd(expr: &Expr, var: &str, ctx: &EvalContext, eps: f64) {
    let sym = expr.diff(var).simplify().eval(ctx).unwrap();
    let mut p = ctx.clone();
    let mut m = ctx.clone();
    let v0 = p.signals.get(var).copied()
        .or_else(|| p.params.get(var).copied())
        .or_else(|| p.behaviors.get(var).copied())
        .unwrap();
    p.signals.insert(var.into(), v0 + eps);
    m.signals.insert(var.into(), v0 - eps);
    if p.params.contains_key(var) {
        p.params.insert(var.into(), v0 + eps);
        m.params.insert(var.into(), v0 - eps);
    }
    let fp = expr.eval(&p).unwrap();
    let fm = expr.eval(&m).unwrap();
    let num = (fp - fm) / (2.0 * eps);
    assert_relative_eq!(sym, num, epsilon = 1e-4);
}

// --- Basic differentiation rules ---
#[test]
fn diff_const_zero() {
    assert_eq!(Expr::Const(5.0).diff("x").simplify(), Expr::Const(0.0));
}

#[test]
fn diff_signal_same() {
    assert_eq!(signal("x").diff("x").simplify(), Expr::Const(1.0));
}

#[test]
fn diff_signal_other() {
    assert_eq!(signal("x").diff("y").simplify(), Expr::Const(0.0));
}

#[test]
fn diff_add() {
    let e = signal("x") + signal("y");
    let d = e.diff("x").simplify();
    assert_eq!(d.eval(&ctx_multi(&[("x", 1.0), ("y", 2.0)])).unwrap(), 1.0);
}

#[test]
fn diff_product_rule() {
    let e = signal("x") * signal("y");
    let ctx = ctx_multi(&[("x", 2.0), ("y", 3.0)]);
    let d = e.diff("x").simplify();
    assert_relative_eq!(d.eval(&ctx).unwrap(), 3.0, epsilon = 1e-9);
}

#[test]
fn diff_product_rule_fd() {
    let e = signal("x") * signal("y");
    let ctx = ctx_multi(&[("x", 1.5), ("y", 2.0)]);
    verify_fd(&e, "x", &ctx, 1e-7);
    verify_fd(&e, "y", &ctx, 1e-7);
}

#[test]
fn diff_x_squared() {
    let x = signal("x");
    let e = x.clone() * x;
    let ctx = ctx_s("x", 3.0);
    let d = e.diff("x").simplify();
    assert_relative_eq!(d.eval(&ctx).unwrap(), 6.0, epsilon = 1e-9);
}

#[test]
fn diff_quotient_fd() {
    let e = signal("x") / signal("y");
    let ctx = ctx_multi(&[("x", 4.0), ("y", 2.0)]);
    verify_fd(&e, "x", &ctx, 1e-7);
    verify_fd(&e, "y", &ctx, 1e-7);
}

#[test]
fn diff_pow_x3() {
    let e = Expr::Pow(Box::new(signal("x")), Box::new(Expr::Const(3.0)));
    let ctx = ctx_s("x", 2.0);
    let d = e.diff("x").simplify();
    assert_relative_eq!(d.eval(&ctx).unwrap(), 12.0, epsilon = 1e-6);
    verify_fd(&e, "x", &ctx, 1e-7);
}

#[test]
fn diff_exp_x() {
    let e = Expr::Exp(Box::new(signal("x")));
    let ctx = ctx_s("x", 1.0);
    let d = e.diff("x").simplify();
    assert_relative_eq!(d.eval(&ctx).unwrap(), std::f64::consts::E, epsilon = 1e-9);
    verify_fd(&e, "x", &ctx, 1e-7);
}

#[test]
fn diff_log_x() {
    let e = Expr::Log(Box::new(signal("x")));
    let ctx = ctx_s("x", 2.0);
    let d = e.diff("x").simplify();
    assert_relative_eq!(d.eval(&ctx).unwrap(), 0.5, epsilon = 1e-9);
    verify_fd(&e, "x", &ctx, 1e-7);
}

#[test]
fn diff_exp_chain_rule() {
    let x = signal("x");
    let e = Expr::Exp(Box::new(x.clone() * x));
    let ctx = ctx_s("x", 0.7);
    verify_fd(&e, "x", &ctx, 1e-7);
}

#[test]
fn diff_log_chain_rule() {
    let x = signal("x");
    let e = Expr::Log(Box::new(Expr::Pow(Box::new(x.clone()), Box::new(Expr::Const(2.0)))));
    let ctx = ctx_s("x", 3.0);
    verify_fd(&e, "x", &ctx, 1e-7);
    let d = e.diff("x").simplify();
    assert_relative_eq!(d.eval(&ctx).unwrap(), 2.0 / 3.0, epsilon = 1e-6);
}

#[test]
fn diff_neg() {
    let e = -signal("x");
    assert_relative_eq!(e.diff("x").simplify().eval(&ctx_s("x", 1.0)).unwrap(), -1.0, epsilon = 1e-12);
}

#[test]
fn diff_constant_base_pow() {
    let e = Expr::Pow(Box::new(Expr::Const(2.0)), Box::new(signal("x")));
    let ctx = ctx_s("x", 3.0);
    verify_fd(&e, "x", &ctx, 1e-7);
}

#[test]
fn diff_both_variable_pow() {
    let e = Expr::Pow(Box::new(signal("x")), Box::new(signal("y")));
    let ctx = ctx_multi(&[("x", 2.0), ("y", 3.0)]);
    verify_fd(&e, "x", &ctx, 1e-7);
    verify_fd(&e, "y", &ctx, 1e-7);
}

// --- Hill derivative ---
#[test]
fn hill_derivative_fd_multiple_points() {
    let e = hill(signal("s"), 5.0, 4.0);
    for s in [1.0, 2.0, 5.0, 10.0, 20.0] {
        verify_fd(&e, "s", &ctx_s("s", s), 1e-7);
    }
}

#[test]
fn hill_derivative_analytical() {
    let s_val = 3.0_f64;
    let h = 5.0_f64;
    let n = 4.0_f64;
    let analytical = n * s_val.powf(n - 1.0) * h.powf(n) / (h.powf(n) + s_val.powf(n)).powi(2);
    let e = hill(signal("s"), h, n);
    let d = e.diff("s").simplify();
    assert_relative_eq!(d.eval(&ctx_s("s", s_val)).unwrap(), analytical, epsilon = 1e-9);
}

// --- Higher-order derivatives ---
#[test]
fn diff_2_x3() {
    let e = Expr::Pow(Box::new(signal("x")), Box::new(Expr::Const(3.0)));
    let d2 = e.diff_n("x", 2);
    let ctx = ctx_s("x", 2.0);
    assert_relative_eq!(d2.eval(&ctx).unwrap(), 12.0, epsilon = 1e-4);
}

#[test]
fn diff_2_exp() {
    let e = Expr::Exp(Box::new(signal("x")));
    let d2 = e.diff_n("x", 2);
    let ctx = ctx_s("x", 1.0);
    assert_relative_eq!(d2.eval(&ctx).unwrap(), std::f64::consts::E, epsilon = 1e-6);
}

#[test]
fn diff_3_x4() {
    let e = Expr::Pow(Box::new(signal("x")), Box::new(Expr::Const(4.0)));
    let d3 = e.diff_n("x", 3);
    let ctx = ctx_s("x", 1.0);
    assert_relative_eq!(d3.eval(&ctx).unwrap(), 24.0, epsilon = 1e-3);
}

// --- Full single-rule derivative ---
#[test]
fn single_rule_derivative_ds() {
    let b0 = 0.0001;
    let bm = 0.0005;
    let h = 5.0;
    let n = 4.0;
    let e = const_(b0) + (const_(bm) - const_(b0)) * hill(signal("s"), h, n);
    let ctx = ctx_s("s", 3.0);
    verify_fd(&e, "s", &ctx, 1e-7);
}

// --- Gradient ---
#[test]
fn gradient_x2_y2() {
    let e = Expr::Pow(Box::new(signal("x")), Box::new(Expr::Const(2.0)))
        + Expr::Pow(Box::new(signal("y")), Box::new(Expr::Const(2.0)));
    let g = e.gradient(&["x", "y"]);
    assert_eq!(g.len(), 2);
    let ctx = ctx_multi(&[("x", 3.0), ("y", 4.0)]);
    assert_relative_eq!(g[0].eval(&ctx).unwrap(), 6.0, epsilon = 1e-6);
    assert_relative_eq!(g[1].eval(&ctx).unwrap(), 8.0, epsilon = 1e-6);
}

// --- Simplification ---
#[test]
fn simplify_zero_plus_x() {
    let x = signal("x");
    assert!(x.structural_eq(&(Expr::Const(0.0) + x.clone()).simplify()));
}

#[test]
fn simplify_x_plus_zero() {
    let x = signal("x");
    assert!(x.structural_eq(&(x.clone() + Expr::Const(0.0)).simplify()));
}

#[test]
fn simplify_x_times_one() {
    let x = signal("x");
    assert!(x.structural_eq(&(x.clone() * Expr::Const(1.0)).simplify()));
}

#[test]
fn simplify_one_times_x() {
    let x = signal("x");
    assert!(x.structural_eq(&(Expr::Const(1.0) * x.clone()).simplify()));
}

#[test]
fn simplify_x_times_zero() {
    assert_eq!((signal("x") * Expr::Const(0.0)).simplify(), Expr::Const(0.0));
}

#[test]
fn simplify_x_sub_x() {
    let x = signal("x");
    assert_eq!((x.clone() - x).simplify(), Expr::Const(0.0));
}

#[test]
fn simplify_x_div_x() {
    let x = signal("x");
    assert_eq!((x.clone() / x).simplify(), Expr::Const(1.0));
}

#[test]
fn simplify_const_fold_add() {
    assert_eq!((Expr::Const(2.0) + Expr::Const(3.0)).simplify(), Expr::Const(5.0));
}

#[test]
fn simplify_const_fold_mul() {
    assert_eq!((Expr::Const(2.0) * Expr::Const(4.0)).simplify(), Expr::Const(8.0));
}

#[test]
fn simplify_double_neg() {
    let x = signal("x");
    assert!(x.structural_eq(&(-(-x.clone())).simplify()));
}

#[test]
fn simplify_log_exp() {
    let x = signal("x");
    let e = Expr::Log(Box::new(Expr::Exp(Box::new(x.clone()))));
    assert!(x.structural_eq(&e.simplify()));
}

#[test]
fn simplify_exp_log() {
    let x = signal("x");
    let e = Expr::Exp(Box::new(Expr::Log(Box::new(x.clone()))));
    assert!(x.structural_eq(&e.simplify()));
}

#[test]
fn simplify_pow_zero() {
    let e = Expr::Pow(Box::new(signal("x")), Box::new(Expr::Const(0.0)));
    assert_eq!(e.simplify(), Expr::Const(1.0));
}

#[test]
fn simplify_pow_one() {
    let x = signal("x");
    let e = Expr::Pow(Box::new(x.clone()), Box::new(Expr::Const(1.0)));
    assert!(x.structural_eq(&e.simplify()));
}

#[test]
fn simplify_neg_const() {
    assert_eq!((-Expr::Const(5.0)).simplify(), Expr::Const(-5.0));
}

#[test]
fn simplify_idempotent() {
    let e = (signal("x") + Expr::Const(0.0)) * Expr::Const(1.0);
    let s1 = e.simplify();
    let s2 = s1.simplify();
    assert!(s1.structural_eq(&s2));
}

#[test]
fn simplify_preserves_eval() {
    let e = (signal("x") * Expr::Const(1.0)) + (Expr::Const(0.0) * signal("y"));
    let ctx = ctx_multi(&[("x", 3.0), ("y", 7.0)]);
    let v1 = e.eval(&ctx).unwrap();
    let v2 = e.simplify().eval(&ctx).unwrap();
    assert_relative_eq!(v1, v2, epsilon = 1e-12);
}

// --- Diff of expression with no var occurrences ---
#[test]
fn diff_no_variable_simplifies_to_zero() {
    let e = Expr::Const(42.0) * Expr::Const(2.0);
    assert_eq!(e.diff("x").simplify(), Expr::Const(0.0));
}

// --- Deep chain rule ---
#[test]
fn diff_deep_chain() {
    let mut e = signal("x");
    for _ in 0..3 {
        e = Expr::Exp(Box::new(e));
    }
    let ctx = ctx_s("x", 0.01);
    verify_fd(&e, "x", &ctx, 1e-7);
}

// --- Stress: variable appears many times ---
#[test]
fn diff_many_variable_occurrences() {
    let x = signal("x");
    let mut e = x.clone();
    for _ in 0..20 {
        e = e + x.clone();
    }
    let ctx = ctx_s("x", 1.0);
    let d = e.diff("x").simplify();
    assert_relative_eq!(d.eval(&ctx).unwrap(), 21.0, epsilon = 1e-6);
}

// --- Product + quotient composite ---
#[test]
fn diff_x2_times_exp_fd() {
    let x = signal("x");
    let e = Expr::Pow(Box::new(x.clone()), Box::new(Expr::Const(2.0)))
        * Expr::Exp(Box::new(x));
    let ctx = ctx_s("x", 1.0);
    verify_fd(&e, "x", &ctx, 1e-7);
}
