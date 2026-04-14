use approx::assert_relative_eq;
use bio_math::expr::{behavior, const_, param, EvalContext, Expr};
use bio_math::integrate::{Euler, Integrator, Rk4};
use bio_math::rule::OdeRule;

fn exp_growth_ode(k: f64) -> OdeRule {
    OdeRule {
        cell_type: "cell".into(),
        behavior: "n".into(),
        rate_expr: Expr::Const(k) * behavior("n"),
        bounds: None,
    }
}

fn ctx_for(beh_name: &str, beh_val: f64) -> EvalContext {
    let mut c = EvalContext::default();
    c.behaviors.insert(beh_name.into(), beh_val);
    c
}

// --- Exponential growth: Euler ---
#[test]
fn euler_exp_growth() {
    let k = 0.1;
    let ode = exp_growth_ode(k);
    let dt = 0.001;
    let steps = 1000;
    let mut y = 1.0;
    for _ in 0..steps {
        y = Euler.step(&ode, &ctx_for("n", y), y, dt).unwrap();
    }
    let t = dt * steps as f64;
    let analytical = (k * t).exp();
    assert_relative_eq!(y, analytical, epsilon = 1e-3);
}

// --- Exponential growth: RK4 ---
#[test]
fn rk4_exp_growth() {
    let k = 0.1;
    let ode = exp_growth_ode(k);
    let dt = 0.01;
    let steps = 100;
    let mut y = 1.0;
    for _ in 0..steps {
        y = Rk4.step(&ode, &ctx_for("n", y), y, dt).unwrap();
    }
    let t = dt * steps as f64;
    let analytical = (k * t).exp();
    assert_relative_eq!(y, analytical, epsilon = 1e-8);
}

// --- Euler convergence O(h) ---
#[test]
fn euler_convergence_order() {
    let k = 0.2;
    let ode = exp_growth_ode(k);
    let t_end = 1.0;
    let analytical = (k * t_end).exp();

    let run = |dt: f64| -> f64 {
        let steps = (t_end / dt).round() as usize;
        let mut y = 1.0;
        for _ in 0..steps {
            y = Euler.step(&ode, &ctx_for("n", y), y, dt).unwrap();
        }
        (y - analytical).abs()
    };
    let e1 = run(0.01);
    let e2 = run(0.005);
    let ratio = e1 / e2;
    assert!(ratio > 1.8 && ratio < 2.3, "ratio={ratio}");
}

// --- RK4 convergence O(h^4) ---
#[test]
fn rk4_convergence_order() {
    let k = 0.2;
    let ode = exp_growth_ode(k);
    let t_end = 1.0;
    let analytical = (k * t_end).exp();

    let run = |dt: f64| -> f64 {
        let steps = (t_end / dt).round() as usize;
        let mut y = 1.0;
        for _ in 0..steps {
            y = Rk4.step(&ode, &ctx_for("n", y), y, dt).unwrap();
        }
        (y - analytical).abs()
    };
    let e1 = run(0.02);
    let e2 = run(0.01);
    let ratio = e1 / e2;
    assert!(ratio > 14.0 && ratio < 18.0, "expected ~16, ratio={ratio}");
}

// --- Logistic growth: RK4 ---
#[test]
fn logistic_rk4() {
    let r = 0.5;
    let cap = 100.0;
    let y0 = 10.0;
    let ode = OdeRule {
        cell_type: "cell".into(),
        behavior: "pop".into(),
        rate_expr: param("r") * behavior("pop") * (const_(1.0) - behavior("pop") / param("K")),
        bounds: None,
    };
    let dt = 0.01;
    let steps = 2000;
    let mut y = y0;
    for _ in 0..steps {
        let mut ctx = EvalContext::default();
        ctx.params.insert("r".into(), r);
        ctx.params.insert("K".into(), cap);
        ctx.behaviors.insert("pop".into(), y);
        y = Rk4.step(&ode, &ctx, y, dt).unwrap();
    }
    let t = dt * steps as f64;
    let analytical = cap / (1.0 + (cap / y0 - 1.0) * (-r * t).exp());
    assert_relative_eq!(y, analytical, epsilon = 1e-4);
}

// --- Bounds clamping ---
#[test]
fn bounds_clamping_upper() {
    let ode = OdeRule {
        cell_type: "cell".into(),
        behavior: "n".into(),
        rate_expr: const_(1000.0),
        bounds: Some((0.0, 50.0)),
    };
    let y = Euler.step(&ode, &ctx_for("n", 40.0), 40.0, 1.0).unwrap();
    assert_eq!(y, 50.0);
}

#[test]
fn bounds_clamping_lower() {
    let ode = OdeRule {
        cell_type: "cell".into(),
        behavior: "n".into(),
        rate_expr: const_(-1000.0),
        bounds: Some((0.0, 50.0)),
    };
    let y = Euler.step(&ode, &ctx_for("n", 10.0), 10.0, 1.0).unwrap();
    assert_eq!(y, 0.0);
}

// --- Single Euler step ---
#[test]
fn euler_single_step() {
    let ode = OdeRule {
        cell_type: "cell".into(),
        behavior: "n".into(),
        rate_expr: const_(2.0),
        bounds: None,
    };
    let y = Euler.step(&ode, &ctx_for("n", 5.0), 5.0, 0.1).unwrap();
    assert_relative_eq!(y, 5.2, epsilon = 1e-12);
}

// --- Single RK4 step ---
#[test]
fn rk4_single_step() {
    let ode = OdeRule {
        cell_type: "cell".into(),
        behavior: "n".into(),
        rate_expr: const_(2.0),
        bounds: None,
    };
    let y = Rk4.step(&ode, &ctx_for("n", 5.0), 5.0, 0.1).unwrap();
    assert_relative_eq!(y, 5.2, epsilon = 1e-12);
}
