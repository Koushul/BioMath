use bio_math::bio_grn;
use bio_math::expr::{behavior, const_, hill, param};
use bio_math::model::CellState;

fn main() {
    let odes = bio_grn! {
        in "cell", d "A" / dt =
            const_(0.2) * hill(param("S"), 0.5, 2.0) * (const_(1.0) - behavior("A"))
            - const_(0.1) * behavior("A"),
            bounded_by (0.0, 1.0);
        in "cell", d "B" / dt =
            const_(0.15) * behavior("A") * (const_(1.0) - behavior("B"))
            - const_(0.08) * behavior("B"),
            bounded_by (0.0, 1.0);
    };
    let m = bio_math::Model::from_rules("two_gene", vec![], Default::default(), None)
        .unwrap()
        .with_ode_rules(odes);

    let mut state = CellState::default();
    state.behaviors.insert("A".into(), 0.1);
    state.behaviors.insert("B".into(), 0.05);
    state.signals.insert("S".into(), 0.8);

    for _ in 0..200 {
        m.step_rk4("cell", &mut state, 0.05).unwrap();
    }
    println!("A={:.3} B={:.3}", state.behaviors["A"], state.behaviors["B"]);
}
