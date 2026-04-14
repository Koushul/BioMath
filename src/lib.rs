pub mod compile;
pub mod dictionary;
pub mod error;
pub mod expr;
pub mod export;
pub mod grammar;
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
pub use grammar::{
    combined_behavior_latex, parse_rule_fields, parse_rule_row, rules_from_csv_str,
    rule_to_grammar_english,
};
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

/// **Model:** optional `rules:` (microenvironment → behavior); optional `paper_csv: ( ... )`
/// (embedded *Cell* 2025-style CSV string); optional `odes:` (coupled ODEs,
/// [`bio_ode_system!`](crate::bio_ode_system) syntax).
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
        $( paper_csv: ( $paper_csv:expr ) )?
        $( odes: { $($ode_body:tt)* } )?
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
        $(
            {
                let __extra = $crate::rules_from_csv_str($paper_csv)
                    .unwrap_or_else(|e| ::core::panic!("paper_csv: {}", e));
                _rules.extend(__extra);
            }
        )?
        let mut _m = $crate::Model::from_rules($name, _rules, _bv, None).expect("model compile");
        $(
            _m = _m.with_ode_rules($crate::bio_ode_system! { $($ode_body)* });
        )?
        _m
    }};
}
