use bio_math::expr::{behavior, const_, hill, param, signal, EvalContext, Expr};
use bio_math::BioMathError;

fn ctx_sig(name: &str, v: f64) -> EvalContext {
    let mut c = EvalContext::default();
    c.signals.insert(name.into(), v);
    c
}

// --- Const ---
#[test]
fn const_zero() {
    assert_eq!(Expr::Const(0.0).eval(&EvalContext::default()).unwrap(), 0.0);
}

#[test]
fn const_negative() {
    assert_eq!(Expr::Const(-3.5).eval(&EvalContext::default()).unwrap(), -3.5);
}

#[test]
fn const_large() {
    assert_eq!(Expr::Const(f64::MAX).eval(&EvalContext::default()).unwrap(), f64::MAX);
}

// --- Add, Sub, Mul, Div, Neg, Pow ---
#[test]
fn add_consts() {
    let e = Expr::Const(2.0) + Expr::Const(3.0);
    assert_eq!(e.eval(&EvalContext::default()).unwrap(), 5.0);
}

#[test]
fn sub_consts() {
    let e = Expr::Const(5.0) - Expr::Const(3.0);
    assert_eq!(e.eval(&EvalContext::default()).unwrap(), 2.0);
}

#[test]
fn mul_consts() {
    let e = Expr::Const(4.0) * Expr::Const(2.5);
    assert_eq!(e.eval(&EvalContext::default()).unwrap(), 10.0);
}

#[test]
fn div_consts() {
    let e = Expr::Const(10.0) / Expr::Const(4.0);
    assert_eq!(e.eval(&EvalContext::default()).unwrap(), 2.5);
}

#[test]
fn div_by_zero() {
    let e = Expr::Const(1.0) / Expr::Const(0.0);
    assert!(matches!(e.eval(&EvalContext::default()), Err(BioMathError::DivisionByZero)));
}

#[test]
fn neg_const() {
    let e = -Expr::Const(7.0);
    assert_eq!(e.eval(&EvalContext::default()).unwrap(), -7.0);
}

#[test]
fn pow_const() {
    let e = Expr::Pow(Box::new(Expr::Const(2.0)), Box::new(Expr::Const(3.0)));
    assert_eq!(e.eval(&EvalContext::default()).unwrap(), 8.0);
}

#[test]
fn pow_zero_zero() {
    let e = Expr::Pow(Box::new(Expr::Const(0.0)), Box::new(Expr::Const(0.0)));
    assert_eq!(e.eval(&EvalContext::default()).unwrap(), 1.0);
}

#[test]
fn pow_negative_base_fractional_exp() {
    let e = Expr::Pow(Box::new(Expr::Const(-1.0)), Box::new(Expr::Const(0.5)));
    assert!(e.eval(&EvalContext::default()).unwrap().is_nan());
}

// --- Signal, Behavior, Param ---
#[test]
fn signal_present() {
    assert_eq!(signal("x").eval(&ctx_sig("x", 3.5)).unwrap(), 3.5);
}

#[test]
fn signal_missing() {
    assert!(matches!(signal("x").eval(&EvalContext::default()), Err(BioMathError::UnknownSignal { .. })));
}

#[test]
fn behavior_present() {
    let mut c = EvalContext::default();
    c.behaviors.insert("b".into(), 42.0);
    assert_eq!(behavior("b").eval(&c).unwrap(), 42.0);
}

#[test]
fn behavior_missing() {
    assert!(matches!(behavior("b").eval(&EvalContext::default()), Err(BioMathError::UnknownBehavior { .. })));
}

#[test]
fn param_present() {
    let mut c = EvalContext::default();
    c.params.insert("k".into(), 0.1);
    assert_eq!(param("k").eval(&c).unwrap(), 0.1);
}

#[test]
fn param_missing() {
    assert!(matches!(param("k").eval(&EvalContext::default()), Err(BioMathError::UnknownParam { .. })));
}

// --- Hill ---
#[test]
fn hill_at_zero() {
    assert_eq!(hill(signal("s"), 5.0, 4.0).eval(&ctx_sig("s", 0.0)).unwrap(), 0.0);
}

#[test]
fn hill_at_half_max() {
    let v = hill(signal("s"), 5.0, 4.0).eval(&ctx_sig("s", 5.0)).unwrap();
    assert!((v - 0.5).abs() < 1e-12);
}

#[test]
fn hill_at_large_s() {
    let v = hill(signal("s"), 5.0, 4.0).eval(&ctx_sig("s", 1e6)).unwrap();
    assert!((v - 1.0).abs() < 1e-9);
}

#[test]
fn hill_various_n() {
    for n in [1.0, 2.0, 4.0, 8.0] {
        let v = hill(signal("s"), 5.0, n).eval(&ctx_sig("s", 5.0)).unwrap();
        assert!((v - 0.5).abs() < 1e-12, "n={n}: expected 0.5, got {v}");
    }
}

// --- Exp, Log ---
#[test]
fn exp_zero() {
    let e = Expr::Exp(Box::new(Expr::Const(0.0)));
    assert!((e.eval(&EvalContext::default()).unwrap() - 1.0).abs() < 1e-15);
}

#[test]
fn exp_one() {
    let e = Expr::Exp(Box::new(Expr::Const(1.0)));
    assert!((e.eval(&EvalContext::default()).unwrap() - std::f64::consts::E).abs() < 1e-12);
}

#[test]
fn log_one() {
    let e = Expr::Log(Box::new(Expr::Const(1.0)));
    assert_eq!(e.eval(&EvalContext::default()).unwrap(), 0.0);
}

#[test]
fn log_e() {
    let e = Expr::Log(Box::new(Expr::Const(std::f64::consts::E)));
    assert!((e.eval(&EvalContext::default()).unwrap() - 1.0).abs() < 1e-12);
}

#[test]
fn log_zero_error() {
    let e = Expr::Log(Box::new(Expr::Const(0.0)));
    assert!(matches!(e.eval(&EvalContext::default()), Err(BioMathError::LogDomain)));
}

#[test]
fn log_negative_error() {
    let e = Expr::Log(Box::new(Expr::Const(-1.0)));
    assert!(matches!(e.eval(&EvalContext::default()), Err(BioMathError::LogDomain)));
}

// --- Abs, Min, Max, Clamp ---
#[test]
fn abs_positive() {
    assert_eq!(Expr::Abs(Box::new(Expr::Const(5.0))).eval(&EvalContext::default()).unwrap(), 5.0);
}

#[test]
fn abs_negative() {
    assert_eq!(Expr::Abs(Box::new(Expr::Const(-5.0))).eval(&EvalContext::default()).unwrap(), 5.0);
}

#[test]
fn abs_zero() {
    assert_eq!(Expr::Abs(Box::new(Expr::Const(0.0))).eval(&EvalContext::default()).unwrap(), 0.0);
}

#[test]
fn min_eval() {
    let e = Expr::Min(Box::new(Expr::Const(3.0)), Box::new(Expr::Const(7.0)));
    assert_eq!(e.eval(&EvalContext::default()).unwrap(), 3.0);
}

#[test]
fn max_eval() {
    let e = Expr::Max(Box::new(Expr::Const(3.0)), Box::new(Expr::Const(7.0)));
    assert_eq!(e.eval(&EvalContext::default()).unwrap(), 7.0);
}

#[test]
fn clamp_below() {
    let e = Expr::Clamp(Box::new(Expr::Const(-1.0)), Box::new(Expr::Const(0.0)), Box::new(Expr::Const(1.0)));
    assert_eq!(e.eval(&EvalContext::default()).unwrap(), 0.0);
}

#[test]
fn clamp_above() {
    let e = Expr::Clamp(Box::new(Expr::Const(5.0)), Box::new(Expr::Const(0.0)), Box::new(Expr::Const(1.0)));
    assert_eq!(e.eval(&EvalContext::default()).unwrap(), 1.0);
}

#[test]
fn clamp_within() {
    let e = Expr::Clamp(Box::new(Expr::Const(0.5)), Box::new(Expr::Const(0.0)), Box::new(Expr::Const(1.0)));
    assert_eq!(e.eval(&EvalContext::default()).unwrap(), 0.5);
}

// --- IfGreater ---
#[test]
fn if_greater_true() {
    let e = Expr::IfGreater(
        Box::new(Expr::Const(5.0)), Box::new(Expr::Const(3.0)),
        Box::new(Expr::Const(1.0)), Box::new(Expr::Const(0.0)),
    );
    assert_eq!(e.eval(&EvalContext::default()).unwrap(), 1.0);
}

#[test]
fn if_greater_false() {
    let e = Expr::IfGreater(
        Box::new(Expr::Const(2.0)), Box::new(Expr::Const(3.0)),
        Box::new(Expr::Const(1.0)), Box::new(Expr::Const(0.0)),
    );
    assert_eq!(e.eval(&EvalContext::default()).unwrap(), 0.0);
}

// --- Nested ---
#[test]
fn deep_nesting() {
    let mut e = signal("x");
    for _ in 0..10 {
        e = e + Expr::Const(1.0);
    }
    assert_eq!(e.eval(&ctx_sig("x", 0.0)).unwrap(), 10.0);
}

// --- Operator overloading ---
#[test]
fn add_tree_shape() {
    let e = signal("x") + const_(1.0);
    assert!(matches!(e, Expr::Add(_, _)));
}

#[test]
fn sub_not_commutative() {
    let ctx = ctx_sig("x", 5.0);
    let a = signal("x") - const_(2.0);
    let b = const_(2.0) - signal("x");
    assert_eq!(a.eval(&ctx).unwrap(), 3.0);
    assert_eq!(b.eval(&ctx).unwrap(), -3.0);
}

#[test]
fn mul_associative() {
    let ctx = ctx_sig("x", 2.0);
    let e = (signal("x") + const_(1.0)) * const_(3.0);
    assert_eq!(e.eval(&ctx).unwrap(), 9.0);
}

// --- free_variables, substitute, is_const, structural_eq ---
#[test]
fn free_variables_collects_all() {
    let e = signal("x") + param("k") * behavior("b");
    let vars = e.free_variables();
    assert!(vars.contains("x"));
    assert!(vars.contains("k"));
    assert!(vars.contains("b"));
    assert_eq!(vars.len(), 3);
}

#[test]
fn substitute_replaces() {
    let e = signal("x") * signal("x");
    let s = e.substitute("x", &Expr::Const(3.0));
    assert_eq!(s.eval(&EvalContext::default()).unwrap(), 9.0);
}

#[test]
fn is_const_true() {
    assert!(Expr::Const(1.0).is_const());
}

#[test]
fn is_const_false() {
    assert!(!signal("x").is_const());
}

#[test]
fn structural_eq_same() {
    let a = signal("x") + const_(1.0);
    let b = signal("x") + const_(1.0);
    assert!(a.structural_eq(&b));
}

#[test]
fn structural_eq_different() {
    let a = signal("x") + const_(1.0);
    let b = signal("x") + const_(2.0);
    assert!(!a.structural_eq(&b));
}

// --- display ---
#[test]
fn display_math_simple() {
    let e = signal("x") + const_(1.0);
    assert_eq!(e.display_math(), "(x + 1)");
}

#[test]
fn display_latex_hill() {
    let e = hill(signal("s"), 5.0, 4.0);
    let l = e.display_latex();
    assert!(l.contains("frac"), "got: {l}");
}
