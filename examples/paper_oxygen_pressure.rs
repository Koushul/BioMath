use bio_math::bio_model;

fn main() {
    let csv = r#"
tumor ; oxygen ; increases ; cycle entry ; 0.0005 ; 5.0 ; 4 ; 0
tumor ; pressure ; decreases ; cycle entry ; 0.0 ; 0.25 ; 4 ; 0
"#;
    let m = bio_model! {
        name: "oxygen_pressure_cycle",
        base_values: {
            ("tumor", "cycle entry") => 0.0001,
        }
        rules: {}
        paper_csv: (csv)
    };
    println!("{}", m.export_summary());
    println!("--- LaTeX (bilinear b, §5.4) ---");
    print!("{}", m.export_grammar_latex());
}
