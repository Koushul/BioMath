use approx::assert_relative_eq;
use bio_math::expr::{EvalContext, Expr};
use bio_math::response::ResponseFnKind;

// --- Hill properties ---
#[test]
fn hill_zero() {
    let h = ResponseFnKind::Hill { half_max: 5.0, hill_power: 4.0 };
    assert_eq!(h.evaluate(0.0).unwrap(), 0.0);
}

#[test]
fn hill_half_max_any_params() {
    for hm in [1.0, 5.0, 10.0, 100.0] {
        for n in [1.0, 2.0, 4.0, 8.0] {
            let h = ResponseFnKind::Hill { half_max: hm, hill_power: n };
            assert_relative_eq!(h.evaluate(hm).unwrap(), 0.5, epsilon = 1e-12);
        }
    }
}

#[test]
fn hill_large_signal_approaches_one() {
    let h = ResponseFnKind::Hill { half_max: 5.0, hill_power: 4.0 };
    assert!(h.evaluate(1e8).unwrap() > 0.999);
}

#[test]
fn hill_bounded_zero_one() {
    let h = ResponseFnKind::Hill { half_max: 5.0, hill_power: 4.0 };
    for s in [0.0, 0.1, 1.0, 5.0, 10.0, 100.0, 1000.0] {
        let v = h.evaluate(s).unwrap();
        assert!((0.0..=1.0).contains(&v), "s={s}, v={v}");
    }
}

#[test]
fn hill_monotone() {
    let h = ResponseFnKind::Hill { half_max: 5.0, hill_power: 4.0 };
    let mut prev = 0.0;
    for s in [0.0, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 50.0] {
        let v = h.evaluate(s).unwrap();
        assert!(v >= prev, "not monotone at s={s}");
        prev = v;
    }
}

#[test]
fn hill_steeper_with_higher_power() {
    let h2 = ResponseFnKind::Hill { half_max: 5.0, hill_power: 2.0 };
    let h8 = ResponseFnKind::Hill { half_max: 5.0, hill_power: 8.0 };
    let below = 4.0;
    let above = 6.0;
    assert!(h8.evaluate(below).unwrap() < h2.evaluate(below).unwrap());
    assert!(h8.evaluate(above).unwrap() > h2.evaluate(above).unwrap());
}

#[test]
fn hill_to_expr_matches_evaluate() {
    let h = ResponseFnKind::Hill { half_max: 5.0, hill_power: 4.0 };
    for s in [0.0, 1.0, 5.0, 10.0, 50.0] {
        let direct = h.evaluate(s).unwrap();
        let expr = h.to_expr("s");
        let mut ctx = EvalContext::default();
        ctx.signals.insert("s".into(), s);
        let via_expr = expr.eval(&ctx).unwrap();
        assert_relative_eq!(direct, via_expr, epsilon = 1e-12);
    }
}

// --- Linear ---
#[test]
fn linear_below() {
    let l = ResponseFnKind::Linear { min_threshold: 1.0, max_threshold: 3.0 };
    assert_eq!(l.evaluate(0.5).unwrap(), 0.0);
}

#[test]
fn linear_above() {
    let l = ResponseFnKind::Linear { min_threshold: 1.0, max_threshold: 3.0 };
    assert_eq!(l.evaluate(4.0).unwrap(), 1.0);
}

#[test]
fn linear_midpoint() {
    let l = ResponseFnKind::Linear { min_threshold: 1.0, max_threshold: 3.0 };
    assert_relative_eq!(l.evaluate(2.0).unwrap(), 0.5, epsilon = 1e-12);
}

#[test]
fn linear_10_points() {
    let l = ResponseFnKind::Linear { min_threshold: 0.0, max_threshold: 10.0 };
    for i in 0..=10 {
        let s = i as f64;
        assert_relative_eq!(l.evaluate(s).unwrap(), s / 10.0, epsilon = 1e-12);
    }
}

// --- Step ---
#[test]
fn step_below() {
    let s = ResponseFnKind::Step { threshold: 5.0 };
    assert_eq!(s.evaluate(4.9).unwrap(), 0.0);
}

#[test]
fn step_at_threshold() {
    let s = ResponseFnKind::Step { threshold: 5.0 };
    assert_eq!(s.evaluate(5.0).unwrap(), 1.0);
}

#[test]
fn step_above() {
    let s = ResponseFnKind::Step { threshold: 5.0 };
    assert_eq!(s.evaluate(5.1).unwrap(), 1.0);
}

// --- Custom ---
#[test]
fn custom_matches_evaluate() {
    let custom_expr = Expr::Hill(
        Box::new(Expr::Signal("s".into())),
        Box::new(Expr::Const(5.0)),
        Box::new(Expr::Const(4.0)),
    );
    let custom = ResponseFnKind::Custom { expr: custom_expr };
    let hill = ResponseFnKind::Hill { half_max: 5.0, hill_power: 4.0 };
    for s in [0.0, 1.0, 5.0, 10.0] {
        assert_relative_eq!(
            custom.evaluate(s).unwrap(),
            hill.evaluate(s).unwrap(),
            epsilon = 1e-12
        );
    }
}
