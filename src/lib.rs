pub mod compile;
pub mod dictionary;
pub mod error;
pub mod expr;
pub mod export;
pub mod germinal_center_gc;
pub mod integrate;
pub mod model;
pub mod response;
pub mod rule;

pub use compile::{
    behavior_rule_set_expr, compile, BaseValueMap, BehaviorRuleSet, CompiledModel, TransitionGraph,
};
pub use dictionary::{
    BehaviorDef, BehaviorKind, Dictionary, SignalCategory, SignalDef, Warning,
};
pub use error::{BioMathError, CompileError};
pub use expr::{behavior, const_, hill, param, signal, EvalContext, Expr};
pub use integrate::{Euler, Integrator, Rk4};
pub use model::{CellState, EvalRequest, EvalResponse, Model, SignalEnv};
pub use response::{r_from_response, t_from_r, ResponseFnKind};
pub use rule::{OdeRule, Response, Rule, RuleBuilder};

#[macro_export]
macro_rules! response_dir {
    (increases) => {
        $crate::Response::Increases
    };
    (decreases) => {
        $crate::Response::Decreases
    };
}

#[macro_export]
macro_rules! bio_rule {
    (
        in $cell:expr, $signal:expr, $dir:tt, $behavior:expr ,
        with $( $param:ident = $val:expr ),+ $(,)?
    ) => {
        $crate::Rule::builder($cell, $signal, $crate::response_dir!($dir), $behavior)
            $( .$param($val) )+
            .build()
    };
}

/// Ordinary differential equation on a **behavior** (rate `d behavior / dt = rhs`).
///
/// Optional `, bounded_by (lo, hi)` clamps the integrated value after each step (e.g. normalized
/// gene activity in `[0, 1]`).
///
/// **Microenvironment vs parameters:** use [`signal`](crate::expr::signal) for cues present in
/// [`CellState::signals`](crate::model::CellState::signals), or [`param`](crate::expr::param) for
/// the same names—[`Model::step`](crate::model::Model::step) copies signals into both `signals` and
/// `params` on the ODE [`EvalContext`](crate::expr::EvalContext).
#[macro_export]
macro_rules! bio_ode {
    (
        in $cell:expr, d $beh:literal / dt = $rhs:expr
        $( , bounded_by ($blo:expr, $bhi:expr) )?
    ) => {{
        let rate = { $rhs };
        #[allow(unused_mut)]
        let mut __bounds: ::core::option::Option<(f64, f64)> = ::core::option::Option::None;
        $(
            __bounds = ::core::option::Option::Some(($blo, $bhi));
        )?
        $crate::OdeRule {
            cell_type: $cell.into(),
            behavior: $beh.into(),
            rate_expr: rate,
            bounds: __bounds,
        }
    }};
}

/// Semicolon-separated list of [`bio_ode!`](crate::bio_ode) equations (gene–gene / protein–protein
/// coupling, mass-action, etc.).
#[macro_export]
macro_rules! bio_ode_system {
    (
        $(
            in $cell:expr, d $beh:literal / dt = $rhs:expr
            $( , bounded_by ($blo:expr, $bhi:expr) )?
        );*
        $(;)?
    ) => {{
        vec![
            $(
                $crate::bio_ode! {
                    in $cell, d $beh / dt = $rhs
                    $( , bounded_by ($blo, $bhi) )?
                }
            ),*
        ]
    }};
}

/// Alias for [`bio_ode_system!`](crate::bio_ode_system): treat **behaviors** as gene (or protein)
/// activity levels coupled by ODEs.
#[macro_export]
macro_rules! bio_grn {
    ($($t:tt)*) => {
        $crate::bio_ode_system! { $($t)* }
    };
}

#[macro_export]
macro_rules! bio_model {
    (
        name: $name:expr,
        $( base_values: { $( ($bk:expr, $bb:expr) => $bval:expr ),+ $(,)? } )?
        rules: {
            $(
                in $cell:expr, $signal:expr, $dir:tt, $behavior:expr ,
                    with $( $param:ident = $pval:expr ),+ $(,)?;
            )*
        }
        $( odes: {
            $(
                in $ocell:expr, d $obeh:literal / dt = $orhs:expr
                $( , bounded_by ($oblo:expr, $obhi:expr) )?
            );*
            $(;)?
        } )?
    ) => {{
        let mut _rules = Vec::new();
        let mut _bv = $crate::BaseValueMap::default();
        $(
            $(
                _bv.insert($bk, $bb, $bval);
            )+
        )?
        $(
            _rules.push($crate::bio_rule! {
                in $cell, $signal, $dir, $behavior ,
                with $( $param = $pval ),+
            });
        )*
        let mut _m = $crate::Model::from_rules($name, _rules, _bv, None).expect("model compile");
        $(
            _m = _m.with_ode_rules($crate::bio_ode_system! {
                $(
                    in $ocell, d $obeh / dt = $orhs
                    $( , bounded_by ($oblo, $obhi) )?
                );*
            });
        )?
        _m
    }};
}
