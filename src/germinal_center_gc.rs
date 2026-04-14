//! Coarse gene–regulatory dynamics for germinal-center B cells:
//! **FOXO1**, **BCL6**, **AICDA**, **CXCR4**, with cues **Tfh_help** and **CXCL12**.
//!
//! Build ODEs with [`germinal_center_gc_odes!`](crate::germinal_center_gc_odes), then attach:
//! `Model::from_rules(...).unwrap().with_ode_rules(germinal_center_gc_odes!(&params))`.

use crate::compile::BaseValueMap;
use crate::model::Model;

/// Default production / coupling strengths (order-of-magnitude; tune for a specific dataset).
#[derive(Clone, Debug)]
pub struct GcGrnParams {
    pub alpha_f: f64,
    pub alpha_b: f64,
    pub alpha_a: f64,
    pub alpha_x: f64,
    pub delta_f: f64,
    pub delta_b: f64,
    pub delta_a: f64,
    pub delta_x: f64,
    pub k_pi3k: f64,
    pub n_hill: f64,
    pub half_bcl6_on_aid: f64,
    pub half_foxo_on_aid: f64,
    pub half_help_on_foxo: f64,
    pub half_cxcl12_on_cxcr4: f64,
}

impl Default for GcGrnParams {
    fn default() -> Self {
        Self {
            alpha_f: 0.35,
            alpha_b: 0.4,
            alpha_a: 0.45,
            alpha_x: 0.5,
            delta_f: 0.25,
            delta_b: 0.22,
            delta_a: 0.28,
            delta_x: 0.24,
            k_pi3k: 0.55,
            n_hill: 2.0,
            half_bcl6_on_aid: 0.35,
            half_foxo_on_aid: 0.35,
            half_help_on_foxo: 0.4,
            half_cxcl12_on_cxcr4: 0.35,
        }
    }
}

/// Germinal-center GRN ODEs. Pass `&GcGrnParams` (e.g. `&GcGrnParams::default()`).
#[macro_export]
macro_rules! germinal_center_gc_odes {
    ($p:expr) => {{
        let __p: &$crate::germinal_center_gc::GcGrnParams = $p;
        $crate::grn! {
            in "gc_B_cell", d "FOXO1" / dt =
                $crate::const_(__p.alpha_f)
                    * $crate::hill($crate::cue!(Tfh_help), __p.half_help_on_foxo, __p.n_hill)
                    * ($crate::const_(1.0) - $crate::const_(__p.k_pi3k) * $crate::cue!(Tfh_help))
                - $crate::const_(__p.delta_f) * $crate::gene!(FOXO1)
                + $crate::const_(0.08) * $crate::gene!(BCL6) * $crate::gene!(FOXO1)
                    * ($crate::const_(1.0) - $crate::gene!(FOXO1)),
                bounded_by (0.0, 1.0);

            in "gc_B_cell", d "BCL6" / dt =
                $crate::const_(__p.alpha_b) * $crate::gene!(FOXO1) * $crate::gene!(BCL6)
                    * ($crate::const_(1.0) - $crate::gene!(BCL6))
                + $crate::const_(0.12) * $crate::gene!(FOXO1) * ($crate::const_(1.0) - $crate::gene!(BCL6))
                - $crate::const_(__p.delta_b) * $crate::gene!(BCL6),
                bounded_by (0.0, 1.0);

            in "gc_B_cell", d "AICDA" / dt =
                $crate::const_(__p.alpha_a)
                    * ($crate::const_(0.55) * $crate::hill($crate::gene!(BCL6), __p.half_bcl6_on_aid, __p.n_hill)
                        + $crate::const_(0.45) * $crate::hill($crate::gene!(FOXO1), __p.half_foxo_on_aid, __p.n_hill))
                    * ($crate::const_(1.0) - $crate::gene!(AICDA))
                - $crate::const_(__p.delta_a) * $crate::gene!(AICDA),
                bounded_by (0.0, 1.0);

            in "gc_B_cell", d "CXCR4" / dt =
                $crate::const_(__p.alpha_x)
                    * $crate::hill($crate::cue!(CXCL12), __p.half_cxcl12_on_cxcr4, __p.n_hill)
                    * $crate::gene!(FOXO1)
                    * ($crate::const_(1.0) - $crate::gene!(CXCR4))
                - $crate::const_(__p.delta_x) * $crate::gene!(CXCR4),
                bounded_by (0.0, 1.0);
        }
    }};
}

pub fn germinal_center_gc_model(params: &GcGrnParams) -> Model {
    Model::from_rules(
        "germinal_center_FOXO1_BCL6_AICDA_CXCR4",
        vec![],
        BaseValueMap::default(),
        None,
    )
    .expect("empty rules compile")
    .with_ode_rules(germinal_center_gc_odes!(params))
}
