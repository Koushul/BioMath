use bio_math::germinal_center_gc::{germinal_center_gc_grn_model, GcGrnParams};
use bio_math::model::CellState;

fn main() {
    let p = GcGrnParams::default();
    let m = germinal_center_gc_grn_model(&p);
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
    }

    println!(
        "FOXO1={:.3} BCL6={:.3} AICDA={:.3} CXCR4={:.3}",
        state.behaviors["FOXO1"],
        state.behaviors["BCL6"],
        state.behaviors["AICDA"],
        state.behaviors["CXCR4"]
    );
}
