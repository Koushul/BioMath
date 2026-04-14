use approx::assert_relative_eq;
use std::collections::HashMap;

use bio_math::compile::{bilinear, BaseValueMap};
use bio_math::expr::{behavior, const_, param, EvalContext};
use bio_math::integrate::Integrator;
use bio_math::model::{CellState, Model};
use bio_math::response::ResponseFnKind;
use bio_math::rule::{OdeRule, Response, Rule};

fn hill_rule(cell: &str, sig: &str, resp: Response, beh: &str, max_r: f64, hm: f64, n: f64) -> Rule {
    Rule {
        cell_type: cell.into(),
        signal: sig.into(),
        response: resp,
        behavior: beh.into(),
        max_response: max_r,
        response_fn: ResponseFnKind::Hill { half_max: hm, hill_power: n },
        applies_to_dead: false,
        condition: None,
    }
}

fn sigs(pairs: &[(&str, f64)]) -> HashMap<String, f64> {
    pairs.iter().map(|(k, v)| ((*k).into(), *v)).collect()
}

fn simple_model() -> Model {
    let rules = vec![
        hill_rule("tumor", "oxygen", Response::Increases, "cycle", 0.0005, 5.0, 4.0),
    ];
    let mut bv = BaseValueMap::default();
    bv.insert("tumor", "cycle", 0.0001);
    Model::from_rules("test", rules, bv, None).unwrap()
}

// --- Single-rule evaluation ---
#[test]
fn up_rule_at_zero_signal() {
    let m = simple_model();
    let b = m.evaluate("tumor", &sigs(&[("oxygen", 0.0)]));
    assert_relative_eq!(b["cycle"], 0.0001, epsilon = 1e-12);
}

#[test]
fn up_rule_at_half_max() {
    let m = simple_model();
    let b = m.evaluate("tumor", &sigs(&[("oxygen", 5.0)]));
    assert_relative_eq!(b["cycle"], (0.0001 + 0.0005) / 2.0, epsilon = 1e-6);
}

#[test]
fn up_rule_at_large_signal() {
    let m = simple_model();
    let b = m.evaluate("tumor", &sigs(&[("oxygen", 1e6)]));
    assert_relative_eq!(b["cycle"], 0.0005, epsilon = 1e-6);
}

#[test]
fn down_rule_at_zero_signal() {
    let rules = vec![
        hill_rule("tumor", "oxygen", Response::Decreases, "necrosis", 0.0, 5.0, 4.0),
    ];
    let mut bv = BaseValueMap::default();
    bv.insert("tumor", "necrosis", 0.01);
    let m = Model::from_rules("test", rules, bv, None).unwrap();
    let b = m.evaluate("tumor", &sigs(&[("oxygen", 0.0)]));
    assert_relative_eq!(b["necrosis"], 0.01, epsilon = 1e-12);
}

#[test]
fn down_rule_at_half_max() {
    let rules = vec![
        hill_rule("tumor", "oxygen", Response::Decreases, "necrosis", 0.0, 5.0, 4.0),
    ];
    let mut bv = BaseValueMap::default();
    bv.insert("tumor", "necrosis", 0.01);
    let m = Model::from_rules("test", rules, bv, None).unwrap();
    let b = m.evaluate("tumor", &sigs(&[("oxygen", 5.0)]));
    assert_relative_eq!(b["necrosis"], (0.01 + 0.0) / 2.0, epsilon = 1e-6);
}

#[test]
fn down_rule_at_large_signal() {
    let rules = vec![
        hill_rule("tumor", "oxygen", Response::Decreases, "necrosis", 0.0, 5.0, 4.0),
    ];
    let mut bv = BaseValueMap::default();
    bv.insert("tumor", "necrosis", 0.01);
    let m = Model::from_rules("test", rules, bv, None).unwrap();
    let b = m.evaluate("tumor", &sigs(&[("oxygen", 1e6)]));
    assert_relative_eq!(b["necrosis"], 0.0, epsilon = 1e-6);
}

// --- Multi-rule bilinear limiting cases ---
#[test]
fn bilinear_u0_d0() {
    assert_relative_eq!(bilinear(1.0, 2.0, 0.5, 0.0, 0.0), 1.0, epsilon = 1e-12);
}

#[test]
fn bilinear_u1_d0() {
    assert_relative_eq!(bilinear(1.0, 2.0, 0.5, 1.0, 0.0), 2.0, epsilon = 1e-12);
}

#[test]
fn bilinear_u0_d1() {
    assert_relative_eq!(bilinear(1.0, 2.0, 0.5, 0.0, 1.0), 0.5, epsilon = 1e-12);
}

#[test]
fn bilinear_u1_d1() {
    assert_relative_eq!(bilinear(1.0, 2.0, 0.5, 1.0, 1.0), 0.5, epsilon = 1e-12);
}

// --- Bounds invariant ---
#[test]
fn bounds_invariant_fuzz() {
    let rules = vec![
        hill_rule("tumor", "oxygen", Response::Increases, "cycle", 0.001, 5.0, 4.0),
        hill_rule("tumor", "pressure", Response::Decreases, "cycle", 0.0, 1.0, 4.0),
    ];
    let mut bv = BaseValueMap::default();
    bv.insert("tumor", "cycle", 0.0001);
    let m = Model::from_rules("test", rules, bv, None).unwrap();
    let set = m.behavior_rule_set("tumor", "cycle").unwrap();
    let b_min = set.min_value;
    let b_max = set.max_value;

    for i in 0..100 {
        let o2 = (i as f64) * 0.5;
        let p = (i as f64) * 0.1;
        let b = m.evaluate("tumor", &sigs(&[("oxygen", o2), ("pressure", p)]));
        let v = b["cycle"];
        assert!(v >= b_min - 1e-12 && v <= b_max + 1e-12,
            "out of bounds: v={v}, min={b_min}, max={b_max}, o2={o2}, p={p}");
    }
}

// --- Modularity ---
#[test]
fn adding_behavior_b2_does_not_change_b1() {
    let rules1 = vec![
        hill_rule("tumor", "oxygen", Response::Increases, "cycle", 0.0005, 5.0, 4.0),
    ];
    let mut bv1 = BaseValueMap::default();
    bv1.insert("tumor", "cycle", 0.0001);
    let m1 = Model::from_rules("m1", rules1, bv1, None).unwrap();

    let rules2 = vec![
        hill_rule("tumor", "oxygen", Response::Increases, "cycle", 0.0005, 5.0, 4.0),
        hill_rule("tumor", "pressure", Response::Decreases, "necrosis", 0.0, 1.0, 4.0),
    ];
    let mut bv2 = BaseValueMap::default();
    bv2.insert("tumor", "cycle", 0.0001);
    bv2.insert("tumor", "necrosis", 0.01);
    let m2 = Model::from_rules("m2", rules2, bv2, None).unwrap();

    let env = sigs(&[("oxygen", 3.0), ("pressure", 0.5)]);
    let b1 = m1.evaluate("tumor", &env);
    let b2 = m2.evaluate("tumor", &env);
    assert_relative_eq!(b1["cycle"], b2["cycle"], epsilon = 1e-12);
}

#[test]
fn different_cell_type_independent() {
    let rules = vec![
        hill_rule("tumor", "oxygen", Response::Increases, "cycle", 0.0005, 5.0, 4.0),
        hill_rule("fibroblast", "ECM", Response::Increases, "secretion", 0.01, 1.0, 2.0),
    ];
    let mut bv = BaseValueMap::default();
    bv.insert("tumor", "cycle", 0.0001);
    bv.insert("fibroblast", "secretion", 0.001);
    let m = Model::from_rules("test", rules, bv, None).unwrap();
    let b = m.evaluate("tumor", &sigs(&[("oxygen", 5.0)]));
    assert!(b.contains_key("cycle"));
    assert!(!b.contains_key("secretion"));
}

// --- Two-signal evaluation ---
#[test]
fn evaluate_two_signals() {
    let rules = vec![
        hill_rule("tumor", "oxygen", Response::Increases, "cycle", 0.0005, 5.0, 4.0),
        hill_rule("tumor", "pressure", Response::Decreases, "cycle", 0.0, 1.0, 4.0),
    ];
    let mut bv = BaseValueMap::default();
    bv.insert("tumor", "cycle", 0.0001);
    let m = Model::from_rules("test", rules, bv, None).unwrap();
    let b = m.evaluate("tumor", &sigs(&[("oxygen", 5.0), ("pressure", 1.0)]));
    assert!(b["cycle"] > 0.0);
    assert!(b["cycle"] < 0.0005);
}

// --- Biphasic peaked response ---
#[test]
fn biphasic_peaked_response() {
    let rules = vec![
        hill_rule("tumor", "ECM", Response::Increases, "speed", 2.0, 3.0, 4.0),
        hill_rule("tumor", "ECM", Response::Decreases, "speed", 0.5, 6.0, 4.0),
    ];
    let mut bv = BaseValueMap::default();
    bv.insert("tumor", "speed", 1.0);
    let m = Model::from_rules("test", rules, bv, None).unwrap();
    let at0 = m.evaluate("tumor", &sigs(&[("ECM", 0.0)]))["speed"];
    let at_mid = m.evaluate("tumor", &sigs(&[("ECM", 4.0)]))["speed"];
    let at_high = m.evaluate("tumor", &sigs(&[("ECM", 100.0)]))["speed"];
    assert_relative_eq!(at0, 1.0, epsilon = 1e-6);
    assert!(at_mid > at0, "peak should be > base");
    assert!(at_high < at_mid, "high signal should be < peak");
}

// --- Sensitivity ---
#[test]
fn sensitivity_expr_nontrivial() {
    let m = simple_model();
    let s = m.sensitivity("tumor", "cycle", "oxygen");
    assert!(s.is_some());
    let expr = s.unwrap();
    assert!(!expr.is_const());
}

#[test]
fn sensitivity_at_half_max_is_maximal() {
    let m = simple_model();
    let s_hm = m.sensitivity_at("tumor", "cycle", "oxygen", &sigs(&[("oxygen", 5.0)]), 1e-6).unwrap();
    let s_below = m.sensitivity_at("tumor", "cycle", "oxygen", &sigs(&[("oxygen", 1.0)]), 1e-6).unwrap();
    let s_above = m.sensitivity_at("tumor", "cycle", "oxygen", &sigs(&[("oxygen", 20.0)]), 1e-6).unwrap();
    assert!(s_hm.abs() > s_below.abs());
    assert!(s_hm.abs() > s_above.abs());
}

// --- ODE stepping ---
#[test]
fn ode_exponential_growth_rk4() {
    let k = 0.1;
    let y0 = 1.0;
    let ode = OdeRule {
        cell_type: "cell".into(),
        behavior: "n".into(),
        rate_expr: param("k") * behavior("n"),
        bounds: None,
    };
    let rules: Vec<Rule> = vec![];
    let m = Model::from_rules("exp", rules, BaseValueMap::default(), None)
        .unwrap()
        .with_ode_rules(vec![ode]);
    let mut state = CellState::default();
    state.behaviors.insert("n".into(), y0);
    state.signals = HashMap::new();
    let mut ctx_params = EvalContext::default();
    ctx_params.params.insert("k".into(), k);

    let dt = 0.001;
    let steps = 1000;
    for _ in 0..steps {
        let cur = state.behaviors["n"];
        let mut ctx = EvalContext::default();
        ctx.params.insert("k".into(), k);
        ctx.behaviors.insert("n".into(), cur);
        let ode_rule = &m.ode_rules[0];
        let next = bio_math::Rk4.step(ode_rule, &ctx, cur, dt).unwrap();
        state.behaviors.insert("n".into(), next);
    }
    let t = dt * steps as f64;
    let analytical = y0 * (k * t).exp();
    assert_relative_eq!(state.behaviors["n"], analytical, epsilon = 1e-6);
}

#[test]
fn ode_logistic_growth_rk4() {
    let r = 0.5;
    let capacity = 100.0;
    let y0 = 10.0;
    let ode = OdeRule {
        cell_type: "cell".into(),
        behavior: "pop".into(),
        rate_expr: param("r") * behavior("pop") * (const_(1.0) - behavior("pop") / param("K")),
        bounds: None,
    };
    let m = Model::from_rules("logistic", vec![], BaseValueMap::default(), None)
        .unwrap()
        .with_ode_rules(vec![ode]);

    let dt = 0.01;
    let steps = 2000;
    let mut cur = y0;
    for _ in 0..steps {
        let mut ctx = EvalContext::default();
        ctx.params.insert("r".into(), r);
        ctx.params.insert("K".into(), capacity);
        ctx.behaviors.insert("pop".into(), cur);
        cur = bio_math::Rk4.step(&m.ode_rules[0], &ctx, cur, dt).unwrap();
    }
    let t = dt * steps as f64;
    let analytical = capacity / (1.0 + (capacity / y0 - 1.0) * (-r * t).exp());
    assert_relative_eq!(cur, analytical, epsilon = 1e-4);
}

// --- Euler convergence order ---
#[test]
fn euler_convergence_order() {
    let k: f64 = 0.1;
    let y0: f64 = 1.0;
    let ode = OdeRule {
        cell_type: "cell".into(),
        behavior: "n".into(),
        rate_expr: param("k") * behavior("n"),
        bounds: None,
    };
    let t_end: f64 = 1.0;
    let analytical = y0 * (k * t_end).exp();

    let run = |dt: f64| -> f64 {
        let steps = (t_end / dt) as usize;
        let mut cur = y0;
        for _ in 0..steps {
            let mut ctx = EvalContext::default();
            ctx.params.insert("k".into(), k);
            ctx.behaviors.insert("n".into(), cur);
            cur = bio_math::Euler.step(&ode, &ctx, cur, dt).unwrap();
        }
        (cur - analytical).abs()
    };

    let e1 = run(0.01);
    let e2 = run(0.005);
    let ratio = e1 / e2;
    assert!(ratio > 1.8 && ratio < 2.2, "Euler O(h): ratio={ratio}");
}

// --- Bounds clamping ---
#[test]
fn ode_bounds_clamping() {
    let ode = OdeRule {
        cell_type: "cell".into(),
        behavior: "n".into(),
        rate_expr: const_(100.0),
        bounds: Some((0.0, 10.0)),
    };
    let mut ctx = EvalContext::default();
    ctx.behaviors.insert("n".into(), 5.0);
    let next = bio_math::Euler.step(&ode, &ctx, 5.0, 1.0).unwrap();
    assert_eq!(next, 10.0);
}

// --- Jacobian ---
#[test]
fn jacobian_2d_system() {
    let ode1 = OdeRule {
        cell_type: "cell".into(),
        behavior: "x".into(),
        rate_expr: behavior("x") - behavior("y"),
        bounds: None,
    };
    let ode2 = OdeRule {
        cell_type: "cell".into(),
        behavior: "y".into(),
        rate_expr: behavior("x") + behavior("y"),
        bounds: None,
    };
    let m = Model::from_rules("j", vec![], BaseValueMap::default(), None)
        .unwrap()
        .with_ode_rules(vec![ode1, ode2]);
    let j = m.jacobian("cell");
    assert_eq!(j.len(), 2);
    assert_eq!(j[0].len(), 2);
    let state = CellState { signals: HashMap::new(), behaviors: sigs(&[("x", 1.0), ("y", 2.0)]) };
    let jv = m.jacobian_at("cell", &state).unwrap();
    assert_relative_eq!(jv[0][0], 1.0, epsilon = 1e-6);
    assert_relative_eq!(jv[0][1], -1.0, epsilon = 1e-6);
    assert_relative_eq!(jv[1][0], 1.0, epsilon = 1e-6);
    assert_relative_eq!(jv[1][1], 1.0, epsilon = 1e-6);
}

// --- Pooled U single rule matches Hill ---
#[test]
fn pooled_u_single_matches_hill() {
    let rules = vec![
        hill_rule("tumor", "oxygen", Response::Increases, "cycle", 0.0005, 5.0, 4.0),
    ];
    let mut bv = BaseValueMap::default();
    bv.insert("tumor", "cycle", 0.0);
    let m = Model::from_rules("test", rules, bv, None).unwrap();
    let b = m.evaluate("tumor", &sigs(&[("oxygen", 5.0)]));
    assert_relative_eq!(b["cycle"], 0.0005 * 0.5, epsilon = 1e-6);
}
