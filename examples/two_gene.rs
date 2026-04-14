use bio_math::model::CellState;
use bio_math::{const_, cue, gene, bio_ode_system, hill, Model};

fn main() {
    let odes = bio_ode_system! {
        in "cell", d "A" / dt =
            const_(0.2) * hill(cue!(S), 0.5, 2.0) * (const_(1.0) - gene!(A))
            - const_(0.1) * gene!(A),
            bounded_by (0.0, 1.0);
        in "cell", d "B" / dt =
            const_(0.15) * gene!(A) * (const_(1.0) - gene!(B))
            - const_(0.08) * gene!(B),
            bounded_by (0.0, 1.0);
    };
    let m = Model::from_rules("two_gene", vec![], Default::default(), None)
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
