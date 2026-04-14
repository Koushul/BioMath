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

#[macro_export]
macro_rules! bio_ode {
    (
        in $cell:expr, d $beh:literal / dt = $rhs:expr
    ) => {{
        let rate = { $rhs };
        $crate::OdeRule {
            cell_type: $cell.into(),
            behavior: $beh.into(),
            rate_expr: rate,
            bounds: None,
        }
    }};
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
        $crate::Model::from_rules($name, _rules, _bv, None).expect("model compile")
    }};
}
