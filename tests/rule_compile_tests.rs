use bio_math::compile::{compile, BaseValueMap};
use bio_math::dictionary::Dictionary;
use bio_math::error::CompileError;
use bio_math::response::ResponseFnKind;
use bio_math::rule::{Response, Rule};

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

// --- b_0 resolution ---
#[test]
fn b0_explicit_override() {
    let rules = vec![hill_rule("tumor", "oxygen", Response::Increases, "cycle", 0.001, 5.0, 4.0)];
    let mut bv = BaseValueMap::default();
    bv.insert("tumor", "cycle", 0.0001);
    let c = compile(rules, &bv, None).unwrap();
    assert_eq!(c.behavior_sets[0].base_value, 0.0001);
}

#[test]
fn b0_from_dictionary() {
    let dict = Dictionary::default_physicell();
    let rules = vec![hill_rule("tumor", "oxygen", Response::Increases, "cycle entry", 0.001, 5.0, 4.0)];
    let bv = BaseValueMap::default();
    let c = compile(rules, &bv, Some(&dict)).unwrap();
    assert_eq!(c.behavior_sets[0].base_value, 0.0001);
}

#[test]
fn b0_zero_when_nothing() {
    let rules = vec![hill_rule("tumor", "oxygen", Response::Increases, "custom_beh", 0.001, 5.0, 4.0)];
    let bv = BaseValueMap::default();
    let c = compile(rules, &bv, None).unwrap();
    assert_eq!(c.behavior_sets[0].base_value, 0.0);
}

// --- Rule grouping ---
#[test]
fn two_behaviors_two_sets() {
    let rules = vec![
        hill_rule("tumor", "oxygen", Response::Increases, "cycle", 0.001, 5.0, 4.0),
        hill_rule("tumor", "pressure", Response::Decreases, "necrosis", 0.0, 1.0, 4.0),
    ];
    let mut bv = BaseValueMap::default();
    bv.insert("tumor", "cycle", 0.0);
    bv.insert("tumor", "necrosis", 0.01);
    let c = compile(rules, &bv, None).unwrap();
    assert_eq!(c.behavior_sets.len(), 2);
}

#[test]
fn same_behavior_one_set_two_up_rules() {
    let rules = vec![
        hill_rule("tumor", "oxygen", Response::Increases, "cycle", 0.001, 5.0, 4.0),
        hill_rule("tumor", "EGF", Response::Increases, "cycle", 0.001, 2.0, 2.0),
    ];
    let mut bv = BaseValueMap::default();
    bv.insert("tumor", "cycle", 0.0);
    let c = compile(rules, &bv, None).unwrap();
    assert_eq!(c.behavior_sets.len(), 1);
    assert_eq!(c.behavior_sets[0].up_rules.len(), 2);
}

#[test]
fn different_cell_types_independent() {
    let rules = vec![
        hill_rule("tumor", "oxygen", Response::Increases, "cycle", 0.001, 5.0, 4.0),
        hill_rule("fibroblast", "ECM", Response::Increases, "secretion", 0.01, 1.0, 2.0),
    ];
    let mut bv = BaseValueMap::default();
    bv.insert("tumor", "cycle", 0.0);
    bv.insert("fibroblast", "secretion", 0.0);
    let c = compile(rules, &bv, None).unwrap();
    assert_eq!(c.behavior_sets.len(), 2);
}

// --- b_M / b_m consistency ---
#[test]
fn inconsistent_max_response_error() {
    let rules = vec![
        hill_rule("tumor", "oxygen", Response::Increases, "cycle", 0.1, 5.0, 4.0),
        hill_rule("tumor", "EGF", Response::Increases, "cycle", 0.2, 2.0, 2.0),
    ];
    let mut bv = BaseValueMap::default();
    bv.insert("tumor", "cycle", 0.01);
    let err = compile(rules, &bv, None).unwrap_err();
    assert!(err.iter().any(|e| matches!(e, CompileError::InconsistentMaxResponse { .. })));
}

#[test]
fn consistent_max_response_ok() {
    let rules = vec![
        hill_rule("tumor", "oxygen", Response::Increases, "cycle", 0.1, 5.0, 4.0),
        hill_rule("tumor", "EGF", Response::Increases, "cycle", 0.1, 2.0, 2.0),
    ];
    let mut bv = BaseValueMap::default();
    bv.insert("tumor", "cycle", 0.01);
    assert!(compile(rules, &bv, None).is_ok());
}

#[test]
fn inconsistent_min_response_error() {
    let rules = vec![
        hill_rule("tumor", "pressure", Response::Decreases, "cycle", 0.0, 1.0, 4.0),
        hill_rule("tumor", "ECM", Response::Decreases, "cycle", 0.001, 2.0, 2.0),
    ];
    let mut bv = BaseValueMap::default();
    bv.insert("tumor", "cycle", 0.01);
    let err = compile(rules, &bv, None).unwrap_err();
    assert!(err.iter().any(|e| matches!(e, CompileError::InconsistentMinResponse { .. })));
}

// --- Biphasic (Section 5.6) ---
#[test]
fn biphasic_same_signal_up_and_down_ok() {
    let rules = vec![
        hill_rule("tumor", "ECM", Response::Increases, "speed", 2.0, 3.0, 4.0),
        hill_rule("tumor", "ECM", Response::Decreases, "speed", 0.5, 6.0, 4.0),
    ];
    let mut bv = BaseValueMap::default();
    bv.insert("tumor", "speed", 1.0);
    let c = compile(rules, &bv, None).unwrap();
    assert_eq!(c.behavior_sets[0].up_rules.len(), 1);
    assert_eq!(c.behavior_sets[0].down_rules.len(), 1);
}

// --- Transition behavior parsing ---
#[test]
fn transition_behavior_graph() {
    let rules = vec![hill_rule("tumor", "oxygen", Response::Decreases, "transition to motile_tumor", 0.001, 5.0, 4.0)];
    let mut bv = BaseValueMap::default();
    bv.insert("tumor", "transition to motile_tumor", 0.002);
    let c = compile(rules, &bv, None).unwrap();
    assert!(c.transitions.is_valid_transition("tumor", "motile_tumor"));
}

// --- Validation errors ---
#[test]
fn dictionary_unknown_signal() {
    let dict = Dictionary::default_physicell();
    let rules = vec![hill_rule("tumor", "xyz_unknown", Response::Increases, "cycle entry", 0.1, 5.0, 4.0)];
    let err = compile(rules, &BaseValueMap::default(), Some(&dict)).unwrap_err();
    assert!(err.iter().any(|e| matches!(e, CompileError::UnknownSignal { .. })));
}

#[test]
fn dictionary_unknown_behavior() {
    let dict = Dictionary::default_physicell();
    let rules = vec![hill_rule("tumor", "oxygen", Response::Increases, "xyz_beh", 0.1, 5.0, 4.0)];
    let err = compile(rules, &BaseValueMap::default(), Some(&dict)).unwrap_err();
    assert!(err.iter().any(|e| matches!(e, CompileError::UnknownBehavior { .. })));
}

#[test]
fn negative_half_max_error() {
    let rules = vec![hill_rule("tumor", "oxygen", Response::Increases, "cycle", 0.1, -1.0, 4.0)];
    let err = compile(rules, &BaseValueMap::default(), None).unwrap_err();
    assert!(err.iter().any(|e| matches!(e, CompileError::NegativeHalfMax { .. })));
}

#[test]
fn zero_half_max_error() {
    let rules = vec![hill_rule("tumor", "oxygen", Response::Increases, "cycle", 0.1, 0.0, 4.0)];
    let err = compile(rules, &BaseValueMap::default(), None).unwrap_err();
    assert!(err.iter().any(|e| matches!(e, CompileError::NonPositiveHalfMax { .. })));
}

#[test]
fn negative_hill_power_error() {
    let rules = vec![hill_rule("tumor", "oxygen", Response::Increases, "cycle", 0.1, 5.0, -2.0)];
    let err = compile(rules, &BaseValueMap::default(), None).unwrap_err();
    assert!(err.iter().any(|e| matches!(e, CompileError::NegativeHillPower { .. })));
}

#[test]
fn compile_ok_single_rule() {
    let rules = vec![hill_rule("tumor", "oxygen", Response::Increases, "cycle", 0.0005, 5.0, 4.0)];
    let mut bv = BaseValueMap::default();
    bv.insert("tumor", "cycle", 0.0001);
    assert!(compile(rules, &bv, None).is_ok());
}
