use std::collections::HashMap;

use bio_math::bio_rule;
use bio_math::compile::BaseValueMap;
use bio_math::Model;

fn main() {
    let mut bv = BaseValueMap::default();
    bv.insert("tumor", "necrosis", 0.01);
    bv.insert("tumor", "cycle entry", 0.0001);
    bv.insert("motile_tumor", "migration speed", 0.0);
    bv.insert("tumor", "transformation to motile_tumor", 0.002);
    bv.insert("motile_tumor", "transformation to tumor", 0.0005);

    let rules = vec![
        bio_rule! {
            in "tumor", "oxygen", decreases, "necrosis" ,
            with half_max = 5.0, hill_power = 4.0, max_response = 0.005
        },
        bio_rule! {
            in "tumor", "oxygen", decreases, "transformation to motile_tumor" ,
            with half_max = 5.0, hill_power = 4.0, max_response = 0.001
        },
        bio_rule! {
            in "motile_tumor", "oxygen", increases, "transformation to tumor" ,
            with half_max = 5.0, hill_power = 4.0, max_response = 0.002
        },
        bio_rule! {
            in "tumor", "oxygen", increases, "cycle entry" ,
            with half_max = 5.0, hill_power = 4.0, max_response = 0.0005
        },
        bio_rule! {
            in "motile_tumor", "oxygen", decreases, "migration speed" ,
            with half_max = 5.0, hill_power = 4.0, max_response = 0.0
        },
    ];

    let model = Model::from_rules("hypoxia", rules, bv, None).expect("compile");
    println!("{}", model.export_summary());

    let mut signals = HashMap::new();
    signals.insert("oxygen".into(), 2.0);
    let out = model.evaluate("tumor", &signals);
    println!("tumor behaviors: {:?}", out);
}
