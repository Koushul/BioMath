use bio_math::germinal_center_gc::{germinal_center_gc_grn_model, GcGrnParams};
use bio_math::model::CellState;

fn steady_state_cxcr4(cxcl12: f64) -> f64 {
    let p = GcGrnParams::default();
    let m = germinal_center_gc_grn_model(&p);
    let mut state = CellState::default();
    state.behaviors.insert("FOXO1".into(), 0.5);
    state.behaviors.insert("BCL6".into(), 0.5);
    state.behaviors.insert("AICDA".into(), 0.2);
    state.behaviors.insert("CXCR4".into(), 0.2);
    state.signals.insert("Tfh_help".into(), 0.5);
    state.signals.insert("CXCL12".into(), cxcl12);
    for _ in 0..800 {
        m.step_rk4("gc_B_cell", &mut state, 0.05).unwrap();
    }
    state.behaviors["CXCR4"]
}

#[test]
fn gc_grn_state_stays_bounded() {
    let p = GcGrnParams::default();
    let m = germinal_center_gc_grn_model(&p);
    let mut state = CellState::default();
    state.behaviors.insert("FOXO1".into(), 0.3);
    state.behaviors.insert("BCL6".into(), 0.3);
    state.behaviors.insert("AICDA".into(), 0.1);
    state.behaviors.insert("CXCR4".into(), 0.3);
    state.signals.insert("Tfh_help".into(), 0.5);
    state.signals.insert("CXCL12".into(), 0.5);
    for _ in 0..200 {
        m.step_rk4("gc_B_cell", &mut state, 0.05).unwrap();
        for v in state.behaviors.values() {
            assert!((0.0..=1.0).contains(v), "out of bounds: {v}");
        }
    }
}

#[test]
fn higher_cxcl12_tends_to_raise_cxcr4() {
    let lo = steady_state_cxcr4(0.1);
    let hi = steady_state_cxcr4(0.95);
    assert!(hi > lo, "expected CXCL12 cue to increase CXCR4 steady level: lo={lo} hi={hi}");
}
