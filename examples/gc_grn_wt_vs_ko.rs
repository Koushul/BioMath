use bio_math::germinal_center_gc::{germinal_center_gc_grn_model, GcGrnParams};
use bio_math::model::CellState;

fn run_to_steady(m: &bio_math::Model, foxo1_ko: bool) -> CellState {
    let mut state = CellState::default();
    state.behaviors.insert("FOXO1".into(), 0.2);
    state.behaviors.insert("BCL6".into(), 0.15);
    state.behaviors.insert("AICDA".into(), 0.05);
    state.behaviors.insert("CXCR4".into(), 0.25);
    state.signals.insert("Tfh_help".into(), 0.6);
    state.signals.insert("CXCL12".into(), 0.7);
    let dt = 0.05;
    for _ in 0..400 {
        m.step_rk4("gc_B_cell", &mut state, dt).expect("step");
        if foxo1_ko {
            state.behaviors.insert("FOXO1".into(), 0.0);
        }
    }
    state
}

fn main() {
    let p = GcGrnParams::default();
    let m = germinal_center_gc_grn_model(&p);
    let wt = run_to_steady(&m, false);
    let ko = run_to_steady(&m, true);
    println!("Same cues: Tfh_help=0.6, CXCL12=0.7; dt=0.05 x 400 steps");
    println!(
        "WT   FOXO1={:.3} BCL6={:.3} AICDA={:.3} CXCR4={:.3}",
        wt.behaviors["FOXO1"],
        wt.behaviors["BCL6"],
        wt.behaviors["AICDA"],
        wt.behaviors["CXCR4"]
    );
    println!(
        "FOXO1 KO (FOXO1 forced to 0 each step): BCL6={:.3} AICDA={:.3} CXCR4={:.3}",
        ko.behaviors["BCL6"],
        ko.behaviors["AICDA"],
        ko.behaviors["CXCR4"]
    );
}
