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

/// Gene (or protein) **level** in an ODE right-hand side: reads [`CellState::behaviors`] under this name.
#[macro_export]
macro_rules! gene {
    ($n:ident) => {
        $crate::behavior(::core::stringify!($n))
    };
    ($n:expr) => {
        $crate::behavior($n)
    };
}

/// Extracellular **cue** (or other input): reads [`CellState::signals`] / `param` under this name
/// (see [`Model::step`](crate::model::Model::step)).
#[macro_export]
macro_rules! cue {
    ($n:ident) => {
        $crate::param(::core::stringify!($n))
    };
    ($n:expr) => {
        $crate::param($n)
    };
}

#[macro_export]
macro_rules! response_dir {
    (increases) => {
        $crate::Response::Increases
    };
    (decreases) => {
        $crate::Response::Decreases
    };
}

/// Hill-style rule: “In `cell_type`, `signal` increases/decreases `behavior` …”.
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

/// Same as [`bio_rule!`](crate::bio_rule); reads as a **signal → behavior** link.
#[macro_export]
macro_rules! signal_rule {
    ($($t:tt)*) => {
        $crate::bio_rule! { $($t)* }
    };
}

/// One ODE: `d name / dt = rhs` on a [`Model`](crate::model::Model) behavior (often a gene level).
///
/// **Bounds (optional tail):** `, bounded_by (lo, hi)` or `, clamp01` (= `[0, 1]`).
///
/// **Inputs:** [`gene!`](crate::gene!) / [`behavior`](crate::expr::behavior) for state;
/// [`cue!`](crate::cue!) / [`param`](crate::expr::param) for microenvironment cues
/// ([`Model::step`](crate::model::Model::step) fills both `signals` and `params` from [`CellState::signals`](crate::model::CellState::signals)).
#[macro_export]
macro_rules! bio_ode {
    (
        in $cell:expr, d $beh:literal / dt = $rhs:expr , bounded_by ($blo:expr, $bhi:expr)
    ) => {{
        let rate = { $rhs };
        $crate::OdeRule {
            cell_type: $cell.into(),
            behavior: $beh.into(),
            rate_expr: rate,
            bounds: ::core::option::Option::Some(($blo, $bhi)),
        }
    }};
    (
        in $cell:expr, d $beh:literal / dt = $rhs:expr , clamp01
    ) => {{
        let rate = { $rhs };
        $crate::OdeRule {
            cell_type: $cell.into(),
            behavior: $beh.into(),
            rate_expr: rate,
            bounds: ::core::option::Option::Some((0.0, 1.0)),
        }
    }};
    (
        in $cell:expr, d $beh:literal / dt = $rhs:expr
    ) => {{
        let rate = { $rhs };
        $crate::OdeRule {
            cell_type: $cell.into(),
            behavior: $beh.into(),
            rate_expr: rate,
            bounds: ::core::option::Option::None,
        }
    }};
}

/// Same as [`bio_ode!`](crate::bio_ode).
#[macro_export]
macro_rules! ode {
    ($($t:tt)*) => {
        $crate::bio_ode! { $($t)* }
    };
}

/// Coupled ODEs (e.g. gene regulatory network): semicolon-separated [`bio_ode!`](crate::bio_ode) clauses.
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

/// Same as [`bio_ode_system!`](crate::bio_ode_system); use for gene–gene (or protein) networks.
#[macro_export]
macro_rules! grn {
    ($($t:tt)*) => {
        $crate::bio_ode_system! { $($t)* }
    };
}

/// Same as [`grn!`](crate::grn) (legacy name).
#[macro_export]
macro_rules! bio_grn {
    ($($t:tt)*) => {
        $crate::grn! { $($t)* }
    };
}

/// **Model:** optional `rules:` (microenvironment → behavior); optional `grn:` and/or `odes:`
/// (same semicolon-separated ODE syntax). If both are present, equations are concatenated in order.
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
        $(
            grn: { $($grn_body:tt)* }
        )?
        $(
            odes: { $($ode_body:tt)* }
        )?
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
        let mut __ode_acc: ::std::vec::Vec<$crate::OdeRule> = ::std::vec::Vec::new();
        $(
            __ode_acc.extend($crate::bio_ode_system! { $($grn_body)* });
        )?
        $(
            __ode_acc.extend($crate::bio_ode_system! { $($ode_body)* });
        )?
        if !__ode_acc.is_empty() {
            _m = _m.with_ode_rules(__ode_acc);
        }
        _m
    }};
}
