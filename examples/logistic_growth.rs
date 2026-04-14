use bio_math::bio_ode;
use bio_math::expr::{behavior, const_, param};
use bio_math::integrate::{Integrator, Rk4};
use bio_math::model::CellState;
use bio_math::Model;

fn main() {
    let ode = bio_ode! {
        in "tumor", d "volume" / dt =
            param("r") * behavior("volume") * (const_(1.0) - behavior("volume") / param("K"))
    };
    let mut state = CellState::default();
    state.behaviors.insert("volume".into(), 10.0);
    let mut ctx = bio_math::EvalContext::default();
    ctx.params.insert("r".into(), 0.2);
    ctx.params.insert("K".into(), 1000.0);
    let m = Model::from_rules("logistic", vec![], Default::default(), None)
        .unwrap()
        .with_ode_rules(vec![ode]);
    for _ in 0..50 {
        let mut ec = ctx.clone();
        for (k, v) in &state.behaviors {
            ec.behaviors.insert(k.clone(), *v);
        }
        let cur = state.behaviors["volume"];
        let next = Rk4.step(m.ode_rules.first().unwrap(), &ec, cur, 0.5).unwrap();
        state.behaviors.insert("volume".into(), next);
    }
    println!("volume after ODE steps: {}", state.behaviors["volume"]);
}
