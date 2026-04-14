use serde::{Deserialize, Serialize};
use thiserror::Error;

#[derive(Debug, Error, Clone, PartialEq, Serialize, Deserialize)]
pub enum BioMathError {
    #[error("unknown signal: {name}")]
    UnknownSignal { name: String },
    #[error("unknown behavior: {name}")]
    UnknownBehavior { name: String },
    #[error("unknown parameter: {name}")]
    UnknownParam { name: String },
    #[error("evaluation error: {0}")]
    EvalError(String),
    #[error("division by zero")]
    DivisionByZero,
    #[error("invalid domain for log")]
    LogDomain,
    #[error("CSV error: {0}")]
    Csv(String),
    #[error("JSON error: {0}")]
    Json(String),
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum CompileError {
    UnknownSignal {
        name: String,
        suggestions: Vec<String>,
    },
    UnknownBehavior {
        name: String,
        suggestions: Vec<String>,
    },
    UnknownCellType {
        name: String,
        suggestions: Vec<String>,
    },
    InconsistentMaxResponse {
        cell_type: String,
        behavior: String,
        values: Vec<f64>,
    },
    InconsistentMinResponse {
        cell_type: String,
        behavior: String,
        values: Vec<f64>,
    },
    InvalidBounds {
        cell_type: String,
        behavior: String,
        base: f64,
        max: f64,
        min: f64,
    },
    NegativeHalfMax {
        rule_index: usize,
    },
    NonPositiveHalfMax {
        rule_index: usize,
    },
    NegativeHillPower {
        rule_index: usize,
    },
    InvalidLinearThresholds {
        rule_index: usize,
    },
    InvalidTransitionTarget {
        from: String,
        to: String,
    },
}

impl From<serde_json::Error> for BioMathError {
    fn from(e: serde_json::Error) -> Self {
        BioMathError::Json(e.to_string())
    }
}
