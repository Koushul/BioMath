use bio_math::{
    bio_model, combined_behavior_latex, parse_rule_row, rules_from_csv_str, rule_to_grammar_english,
    Model,
};
use bio_math::compile::{compile, pooled_u_expr, BaseValueMap};

#[test]
fn parse_eight_column_hill() {
    let r = parse_rule_row("tumor ; oxygen ; increases ; cycle ; 0.0005 ; 5.0 ; 4 ; 0").unwrap();
    assert_eq!(r.cell_type, "tumor");
    assert_eq!(r.signal, "oxygen");
    assert!(r.to_english().contains("half-max 5"));
}

#[test]
fn parse_eleven_column_linear() {
    let line = "tumor ; ECM ; increases ; migration speed ; 10.0 ; 0 ; 1 ; 0 ; linear ; 2.0 ; 7.0";
    let r = parse_rule_row(line).unwrap();
    assert!(rule_to_grammar_english(&r).contains("linear"));
}

#[test]
fn csv_roundtrip_row_shape() {
    let r = parse_rule_row("a ; s ; decreases ; necrosis ; 0 ; 5 ; 2 ; 1").unwrap();
    let row = r.to_csv_row();
    let parts: Vec<&str> = row.split(';').map(|x| x.trim()).collect();
    assert_eq!(parts.len(), 11);
}

#[test]
fn paper_csv_in_bio_model() {
    let csv = r#"
# oxygen + pressure on cycle (Cell 2025 style)
tumor ; oxygen ; increases ; cycle ; 0.0005 ; 5.0 ; 4 ; 0
tumor ; pressure ; decreases ; cycle ; 0.0 ; 0.25 ; 4 ; 0
"#;
    let m = bio_model! {
        name: "from_paper_csv",
        base_values: {
            ("tumor", "cycle") => 0.0001,
        }
        rules: {}
        paper_csv: (csv)
    };
    let env = [("oxygen".into(), 5.0), ("pressure".into(), 0.25)].into_iter().collect();
    let b = m.evaluate("tumor", &env);
    assert!(b["cycle"] > 0.0);
    assert!(b["cycle"] < 0.0005);
}

#[test]
fn combined_latex_nonempty() {
    let rules = rules_from_csv_str(
        "tumor ; oxygen ; increases ; cycle ; 0.0005 ; 5.0 ; 4 ; 0\n\
         tumor ; pressure ; decreases ; cycle ; 0.0 ; 0.25 ; 4 ; 0\n",
    )
    .unwrap();
    let compiled = compile(rules, &BaseValueMap::default(), None).unwrap();
    let set = compiled
        .behavior_sets
        .iter()
        .find(|s| s.behavior == "cycle")
        .unwrap();
    let tex = combined_behavior_latex(set);
    assert!(tex.contains("frac"));
}

#[test]
fn pooled_u_latex_single_hill() {
    let rules = rules_from_csv_str("tumor ; oxygen ; increases ; cycle ; 0.0005 ; 5.0 ; 4 ; 0\n")
        .unwrap();
    let compiled = compile(rules, &BaseValueMap::default(), None).unwrap();
    let set = &compiled.behavior_sets[0];
    let u = pooled_u_expr(&set.up_rules).simplify().display_latex();
    assert!(u.contains("frac"));
}

#[test]
fn model_from_csv_linear_row() {
    let csv = "c ; x ; increases ; b ; 1 ; 0 ; 1 ; 0 ; linear ; 0 ; 10\n";
    let m = Model::from_csv_reader(csv.as_bytes()).unwrap();
    assert_eq!(m.compiled.raw_rules.len(), 1);
}
