//! Coarse gene–regulatory dynamics for germinal-center B cells:
//! **FOXO1**, **BCL6**, **AICDA**, **CXCR4**, with external cues **Tfh_help**
//! (CD40 / selection-associated T cell help) and **CXCL12** (SDF-1, dark-zone chemokine).
//!
//! Literature anchors (mechanisms encoded qualitatively, not fitted to time courses):
//! - **FOXO1** promotes GC progression, proliferation, and LZ→DZ transition (e.g. BATF axis);
//!   PI3K–AKT opposes nuclear FOXO1 (*Immunity* 2015; *Nat. Immunol.* / PMC work on GC FOXO1).
//! - **BCL6** is the master GC transcriptional repressor; supports **AID** / *AICDA* in part by
//!   repressing *miR-155* (*Blood* 2012; PMC3526356).
//! - **FOXO1** is required for normal peripheral B cell programs relevant to CSR / **AID**
//!   (*Nat. Immunol.* Foxo1 stages of B cell differentiation).
//! - **CXCR4** enforces dark-zone positioning and DZ/LZ organization with CXCL12 gradients
//!   (*Nat. Immunol.* 2000; subsequent GC B cell work).
//!
//! State variables are dimensionless activity levels in \([0, 1]\) (not absolute mRNA counts).
//! The ODE list is emitted with [`crate::bio_grn!`] so the same macro surface used for other
//! gene–gene models applies here; strengths come from [`GcGrnParams`].

use crate::compile::BaseValueMap;
use crate::expr::{behavior, const_, hill, param};
use crate::model::Model;
use crate::rule::OdeRule;

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

macro_rules! germinal_center_gc_grn_odes_for {
    ($p:ident) => {
        $crate::bio_grn! {
            in "gc_B_cell", d "FOXO1" / dt =
                const_($p.alpha_f) * hill(param("Tfh_help"), $p.half_help_on_foxo, $p.n_hill)
                    * (const_(1.0) - const_($p.k_pi3k) * param("Tfh_help"))
                - const_($p.delta_f) * behavior("FOXO1")
                + const_(0.08) * behavior("BCL6") * behavior("FOXO1")
                    * (const_(1.0) - behavior("FOXO1")),
                bounded_by (0.0, 1.0);

            in "gc_B_cell", d "BCL6" / dt =
                const_($p.alpha_b) * behavior("FOXO1") * behavior("BCL6")
                    * (const_(1.0) - behavior("BCL6"))
                + const_(0.12) * behavior("FOXO1") * (const_(1.0) - behavior("BCL6"))
                - const_($p.delta_b) * behavior("BCL6"),
                bounded_by (0.0, 1.0);

            in "gc_B_cell", d "AICDA" / dt =
                const_($p.alpha_a)
                    * (const_(0.55) * hill(behavior("BCL6"), $p.half_bcl6_on_aid, $p.n_hill)
                        + const_(0.45) * hill(behavior("FOXO1"), $p.half_foxo_on_aid, $p.n_hill))
                    * (const_(1.0) - behavior("AICDA"))
                - const_($p.delta_a) * behavior("AICDA"),
                bounded_by (0.0, 1.0);

            in "gc_B_cell", d "CXCR4" / dt =
                const_($p.alpha_x)
                    * hill(param("CXCL12"), $p.half_cxcl12_on_cxcr4, $p.n_hill)
                    * behavior("FOXO1")
                    * (const_(1.0) - behavior("CXCR4"))
                - const_($p.delta_x) * behavior("CXCR4"),
                bounded_by (0.0, 1.0);
        }
    };
}

pub fn germinal_center_gc_grn_odes(p: &GcGrnParams) -> Vec<OdeRule> {
    germinal_center_gc_grn_odes_for!(p)
}

pub fn germinal_center_gc_grn_model(p: &GcGrnParams) -> Model {
    Model::from_rules(
        "germinal_center_FOXO1_BCL6_AICDA_CXCR4",
        vec![],
        BaseValueMap::default(),
        None,
    )
    .expect("empty rules compile")
    .with_ode_rules(germinal_center_gc_grn_odes(p))
}
