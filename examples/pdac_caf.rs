use std::collections::HashMap;

use bio_math::bio_rule;
use bio_math::compile::BaseValueMap;
use bio_math::Model;

fn main() {
    let mut bv = BaseValueMap::default();
    bv.insert("fibroblast", "ECM secretion", 0.0);
    bv.insert("epithelial_tumor", "cycle entry", 0.0001);
    bv.insert("mesenchymal_tumor", "migration speed", 1.0);

    let rules = vec![
        bio_rule! {
            in "fibroblast", "ECM", increases, "ECM secretion" ,
            with half_max = 1.0, hill_power = 4.0, max_response = 0.01
        },
        bio_rule! {
            in "epithelial_tumor", "fibroblast_factor", increases, "transformation to mesenchymal_tumor" ,
            with half_max = 0.5, hill_power = 2.0, max_response = 0.002
        },
        bio_rule! {
            in "mesenchymal_tumor", "ECM", increases, "migration speed" ,
            with half_max = 3.0, hill_power = 4.0, max_response = 2.0
        },
        bio_rule! {
            in "mesenchymal_tumor", "ECM", decreases, "migration speed" ,
            with half_max = 6.0, hill_power = 4.0, max_response = 0.5
        },
        bio_rule! {
            in "epithelial_tumor", "inflammatory_factor", increases, "cycle entry" ,
            with half_max = 0.25, hill_power = 2.0, max_response = 0.0008
        },
    ];

    let model = Model::from_rules("pdac_caf", rules, bv, None).expect("compile");
    println!("{}", model.export_summary());

    let mut s = HashMap::new();
    s.insert("ECM".into(), 4.0);
    s.insert("fibroblast_factor".into(), 0.3);
    s.insert("inflammatory_factor".into(), 0.2);
    let v = model.evaluate("mesenchymal_tumor", &s);
    println!("mesenchymal_tumor migration speed (biphasic): {:?}", v.get("migration speed"));
}
