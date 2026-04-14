use approx::assert_relative_eq;
use std::collections::HashMap;

use bio_math::{bio_grn, bio_model, bio_ode, bio_ode_system, bio_rule};
use bio_math::expr::{behavior, const_, param, EvalContext};
use bio_math::response::ResponseFnKind;
use bio_math::rule::Response;

// --- bio_rule! increases ---
#[test]
fn bio_rule_increases() {
    let r = bio_rule! {
        in "tumor", "oxygen", increases, "cycle entry" ,
        with half_max = 5.0, hill_power = 4.0, max_response = 0.0005
    };
    assert_eq!(r.cell_type, "tumor");
    assert_eq!(r.signal, "oxygen");
    assert_eq!(r.response, Response::Increases);
    assert_eq!(r.behavior, "cycle entry");
    assert_eq!(r.max_response, 0.0005);
    assert!(matches!(r.response_fn, ResponseFnKind::Hill { half_max, hill_power }
        if (half_max - 5.0).abs() < 1e-12 && (hill_power - 4.0).abs() < 1e-12));
}

// --- bio_rule! decreases ---
#[test]
fn bio_rule_decreases() {
    let r = bio_rule! {
        in "tumor", "pressure", decreases, "cycle entry" ,
        with half_max = 0.25, hill_power = 4.0, max_response = 0.0
    };
    assert_eq!(r.response, Response::Decreases);
    assert_eq!(r.behavior, "cycle entry");
}

// --- bio_rule! field verification ---
#[test]
fn bio_rule_fields() {
    let r = bio_rule! {
        in "fibroblast", "ECM", increases, "migration speed" ,
        with half_max = 1.0, hill_power = 2.0, max_response = 10.0
    };
    assert_eq!(r.cell_type, "fibroblast");
    assert_eq!(r.signal, "ECM");
    assert!(!r.applies_to_dead);
    assert!(r.condition.is_none());
}

// --- bio_ode! basic ---
#[test]
fn bio_ode_compiles_and_evaluates() {
    let ode = bio_ode! {
        in "tumor", d "volume" / dt =
            param("r") * behavior("volume") * (const_(1.0) - behavior("volume") / param("K"))
    };
    assert_eq!(ode.cell_type, "tumor");
    assert_eq!(ode.behavior, "volume");
    let mut ctx = EvalContext::default();
    ctx.params.insert("r".into(), 0.1);
    ctx.params.insert("K".into(), 100.0);
    ctx.behaviors.insert("volume".into(), 50.0);
    let rate = ode.rate_expr.eval(&ctx).unwrap();
    assert_relative_eq!(rate, 0.1 * 50.0 * (1.0 - 50.0 / 100.0), epsilon = 1e-12);
}

// --- bio_ode! with signal ---
#[test]
fn bio_ode_with_signal() {
    let ode = bio_ode! {
        in "cell", d "n" / dt = const_(0.1) * behavior("n")
    };
    assert!(ode.bounds.is_none());
    let mut ctx = EvalContext::default();
    ctx.behaviors.insert("n".into(), 10.0);
    let rate = ode.rate_expr.eval(&ctx).unwrap();
    assert_relative_eq!(rate, 1.0, epsilon = 1e-12);
}

#[test]
fn bio_ode_bounded_by() {
    let ode = bio_ode! {
        in "cell", d "g" / dt = const_(1.0),
        bounded_by (0.0, 0.5)
    };
    assert_eq!(ode.bounds, Some((0.0, 0.5)));
}

#[test]
fn bio_ode_system_two_genes() {
    let odes = bio_ode_system! {
        in "x", d "A" / dt = const_(0.1) * behavior("B") - const_(0.05) * behavior("A"), bounded_by (0.0, 1.0);
        in "x", d "B" / dt = const_(0.1) * behavior("A") - const_(0.05) * behavior("B"), bounded_by (0.0, 1.0);
    };
    assert_eq!(odes.len(), 2);
    assert_eq!(odes[0].behavior, "A");
    assert_eq!(odes[1].behavior, "B");
}

#[test]
fn bio_grn_alias_matches_system() {
    let a = bio_grn! {
        in "c", d "u" / dt = const_(0.0), bounded_by (0.0, 1.0);
    };
    let b = bio_ode_system! {
        in "c", d "u" / dt = const_(0.0), bounded_by (0.0, 1.0);
    };
    assert_eq!(a.len(), b.len());
    assert_eq!(a[0].behavior, b[0].behavior);
}

#[test]
fn bio_model_with_odes_block() {
    let m = bio_model! {
        name: "rules_plus_grn",
        base_values: {
            ("tumor", "cycle") => 0.0001,
        }
        rules: {
            in "tumor", "oxygen", increases, "cycle" ,
                with half_max = 5.0, hill_power = 4.0, max_response = 0.0005;
        }
        odes: {
            in "tumor", d "gene_x" / dt = const_(0.0), bounded_by (0.0, 1.0);
        }
    };
    assert_eq!(m.ode_rules.len(), 1);
    assert_eq!(m.ode_rules[0].behavior, "gene_x");
    let env: HashMap<String, f64> = [("oxygen".into(), 5.0)].into_iter().collect();
    let b = m.evaluate("tumor", &env);
    assert!(b.contains_key("cycle"));
}

// --- bio_model! basic ---
#[test]
fn bio_model_builds_and_evaluates() {
    let model = bio_model! {
        name: "hypoxia",
        base_values: {
            ("tumor", "necrosis") => 0.01,
            ("tumor", "cycle entry") => 0.0001,
        }
        rules: {
            in "tumor", "oxygen", decreases, "necrosis" ,
                with half_max = 5.0, hill_power = 4.0, max_response = 0.005;
            in "tumor", "oxygen", increases, "cycle entry" ,
                with half_max = 5.0, hill_power = 4.0, max_response = 0.0005;
        }
    };
    assert_eq!(model.name, "hypoxia");
    let env: HashMap<String, f64> = [("oxygen".into(), 5.0)].into_iter().collect();
    let b = model.evaluate("tumor", &env);
    assert!(b.contains_key("necrosis"));
    assert!(b.contains_key("cycle entry"));
}

// --- bio_model! hand-computed values ---
#[test]
fn bio_model_at_zero_signal() {
    let model = bio_model! {
        name: "test",
        base_values: {
            ("tumor", "cycle") => 0.0001,
        }
        rules: {
            in "tumor", "oxygen", increases, "cycle" ,
                with half_max = 5.0, hill_power = 4.0, max_response = 0.0005;
        }
    };
    let env: HashMap<String, f64> = [("oxygen".into(), 0.0)].into_iter().collect();
    let b = model.evaluate("tumor", &env);
    assert_relative_eq!(b["cycle"], 0.0001, epsilon = 1e-12);
}

// --- bio_model! multiple rules ---
#[test]
fn bio_model_multiple_rules_semicolons() {
    let model = bio_model! {
        name: "multi",
        base_values: {
            ("tumor", "cycle") => 0.0001,
        }
        rules: {
            in "tumor", "oxygen", increases, "cycle" ,
                with half_max = 5.0, hill_power = 4.0, max_response = 0.0005;
            in "tumor", "pressure", decreases, "cycle" ,
                with half_max = 0.25, hill_power = 4.0, max_response = 0.0;
        }
    };
    let env: HashMap<String, f64> = [("oxygen".into(), 5.0), ("pressure".into(), 0.25)].into_iter().collect();
    let b = model.evaluate("tumor", &env);
    assert!(b["cycle"] > 0.0);
    assert!(b["cycle"] < 0.0005);
}
