use approx::assert_relative_eq;
use std::collections::HashMap;

use bio_math::compile::BaseValueMap;
use bio_math::expr::{behavior, const_, hill, param, signal, EvalContext, Expr};
use bio_math::model::Model;
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

// --- Expr variants roundtrip ---
#[test]
fn expr_const_roundtrip() {
    let e = Expr::Const(42.0);
    let j = serde_json::to_string(&e).unwrap();
    let e2: Expr = serde_json::from_str(&j).unwrap();
    assert_eq!(e, e2);
}

#[test]
fn expr_signal_roundtrip() {
    let e = Expr::Signal("oxygen".into());
    let j = serde_json::to_string(&e).unwrap();
    let e2: Expr = serde_json::from_str(&j).unwrap();
    assert_eq!(e, e2);
}

#[test]
fn expr_add_roundtrip() {
    let e = signal("x") + const_(1.0);
    let j = serde_json::to_string(&e).unwrap();
    let e2: Expr = serde_json::from_str(&j).unwrap();
    assert_eq!(e, e2);
}

#[test]
fn expr_hill_roundtrip() {
    let e = hill(signal("s"), 5.0, 4.0);
    let j = serde_json::to_string(&e).unwrap();
    let e2: Expr = serde_json::from_str(&j).unwrap();
    assert_eq!(e, e2);
}

#[test]
fn expr_nested_roundtrip() {
    let e = Expr::Exp(Box::new(
        signal("x") * signal("x") + Expr::Pow(Box::new(signal("y")), Box::new(const_(2.0)))
    ));
    let j = serde_json::to_string(&e).unwrap();
    let e2: Expr = serde_json::from_str(&j).unwrap();
    assert_eq!(e, e2);
}

#[test]
fn expr_ifgreater_roundtrip() {
    let e = Expr::IfGreater(
        Box::new(signal("x")), Box::new(const_(0.0)),
        Box::new(const_(1.0)), Box::new(const_(-1.0)),
    );
    let j = serde_json::to_string(&e).unwrap();
    let e2: Expr = serde_json::from_str(&j).unwrap();
    assert_eq!(e, e2);
}

#[test]
fn expr_abs_min_max_clamp_roundtrip() {
    for e in [
        Expr::Abs(Box::new(signal("x"))),
        Expr::Min(Box::new(signal("x")), Box::new(const_(5.0))),
        Expr::Max(Box::new(signal("x")), Box::new(const_(5.0))),
        Expr::Clamp(Box::new(signal("x")), Box::new(const_(0.0)), Box::new(const_(10.0))),
    ] {
        let j = serde_json::to_string(&e).unwrap();
        let e2: Expr = serde_json::from_str(&j).unwrap();
        assert_eq!(e, e2);
    }
}

// --- Diff result roundtrip ---
#[test]
fn diff_result_roundtrip_eval() {
    let e = hill(signal("s"), 5.0, 4.0);
    let d = e.diff("s").simplify();
    let j = serde_json::to_string(&d).unwrap();
    let d2: Expr = serde_json::from_str(&j).unwrap();
    let ctx = EvalContext {
        signals: [("s".into(), 3.0)].into_iter().collect(),
        ..Default::default()
    };
    let v1 = d.eval(&ctx).unwrap();
    let v2 = d2.eval(&ctx).unwrap();
    assert_relative_eq!(v1, v2, epsilon = 1e-12);
}

// --- Rule roundtrip ---
#[test]
fn rule_json_roundtrip() {
    let r = hill_rule("tumor", "oxygen", Response::Increases, "cycle", 0.0005, 5.0, 4.0);
    let j = serde_json::to_string(&r).unwrap();
    let r2: Rule = serde_json::from_str(&j).unwrap();
    assert_eq!(r.to_english(), r2.to_english());
    assert_eq!(r, r2);
}

// --- OdeRule roundtrip ---
#[test]
fn ode_rule_json_roundtrip() {
    let ode = OdeRule {
        cell_type: "cell".into(),
        behavior: "volume".into(),
        rate_expr: param("r") * behavior("volume") * (const_(1.0) - behavior("volume") / param("K")),
        bounds: Some((0.0, 1000.0)),
    };
    let j = serde_json::to_string(&ode).unwrap();
    let ode2: OdeRule = serde_json::from_str(&j).unwrap();
    let ctx = EvalContext {
        params: [("r".into(), 0.1), ("K".into(), 100.0)].into_iter().collect(),
        behaviors: [("volume".into(), 50.0)].into_iter().collect(),
        ..Default::default()
    };
    assert_relative_eq!(
        ode.rate_expr.eval(&ctx).unwrap(),
        ode2.rate_expr.eval(&ctx).unwrap(),
        epsilon = 1e-12
    );
}

// --- Model roundtrip ---
#[test]
fn model_json_roundtrip() {
    let rules = vec![
        hill_rule("tumor", "oxygen", Response::Increases, "cycle", 0.0005, 5.0, 4.0),
        hill_rule("tumor", "pressure", Response::Decreases, "cycle", 0.0, 1.0, 4.0),
    ];
    let mut bv = BaseValueMap::default();
    bv.insert("tumor", "cycle", 0.0001);
    let m = Model::from_rules("test", rules, bv, None).unwrap();

    let j = m.to_json().unwrap();
    let m2 = Model::from_json(&j).unwrap();

    for o2 in [0.0, 2.0, 5.0, 10.0, 50.0] {
        for p in [0.0, 0.5, 1.0, 5.0] {
            let env = sigs(&[("oxygen", o2), ("pressure", p)]);
            let b1 = m.evaluate("tumor", &env);
            let b2 = m2.evaluate("tumor", &env);
            assert_relative_eq!(b1["cycle"], b2["cycle"], epsilon = 1e-12);
        }
    }
}

#[test]
fn model_json_roundtrip_with_base_values() {
    let rules = vec![
        hill_rule("tumor", "oxygen", Response::Increases, "cycle", 0.001, 5.0, 4.0),
    ];
    let mut bv = BaseValueMap::default();
    bv.insert("tumor", "cycle", 0.0002);
    let m = Model::from_rules("base_test", rules, bv, None).unwrap();
    let j = m.to_json().unwrap();
    let m2 = Model::from_json(&j).unwrap();
    let env = sigs(&[("oxygen", 5.0)]);
    assert_relative_eq!(m.evaluate("tumor", &env)["cycle"], m2.evaluate("tumor", &env)["cycle"], epsilon = 1e-12);
}

// --- CSV roundtrip ---
#[test]
fn csv_roundtrip() {
    let rules = vec![
        hill_rule("tumor", "oxygen", Response::Increases, "cycle", 0.0005, 5.0, 4.0),
        hill_rule("tumor", "pressure", Response::Decreases, "cycle", 0.0, 1.0, 4.0),
    ];
    let mut bv = BaseValueMap::default();
    bv.insert("tumor", "cycle", 0.0001);
    let m = Model::from_rules("csv_test", rules, bv, None).unwrap();

    let mut csv_buf = Vec::new();
    m.to_csv_writer(&mut csv_buf).unwrap();
    let csv_str = String::from_utf8(csv_buf.clone()).unwrap();
    assert!(csv_str.contains("tumor"));

    let m2 = Model::from_csv_reader(csv_buf.as_slice()).unwrap();
    assert_eq!(m2.compiled.raw_rules.len(), m.compiled.raw_rules.len());
}

// --- EvalRequest/Response ---
#[test]
fn eval_request_json() {
    let req = bio_math::EvalRequest {
        cell_type: "tumor".into(),
        signals: sigs(&[("oxygen", 5.0)]),
        state: HashMap::new(),
        dt: Some(0.01),
    };
    let j = serde_json::to_string(&req).unwrap();
    let req2: bio_math::EvalRequest = serde_json::from_str(&j).unwrap();
    assert_eq!(req2.cell_type, "tumor");
    assert_relative_eq!(req2.signals["oxygen"], 5.0, epsilon = 1e-12);
}
