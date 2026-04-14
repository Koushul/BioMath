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

use crate::expr::{behavior, const_, hill, param};
use crate::model::Model;
use crate::rule::OdeRule;
use crate::compile::BaseValueMap;

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

/// \(\Phi(u) = u^n / (h^n + u^n)\) with `u` a behavior level.
pub fn germinal_center_gc_grn_odes(p: &GcGrnParams) -> Vec<OdeRule> {
    let tfh = param("Tfh_help");
    let cxcl12 = param("CXCL12");
    let foxo = behavior("FOXO1");
    let bcl6 = behavior("BCL6");
    let aicda = behavior("AICDA");
    let cxcr4 = behavior("CXCR4");

    let help_on_foxo = hill(tfh.clone(), p.half_help_on_foxo, p.n_hill);
    let cxcl12_on_cxcr4 = hill(cxcl12.clone(), p.half_cxcl12_on_cxcr4, p.n_hill);
    let bcl6_on_aid = hill(bcl6.clone(), p.half_bcl6_on_aid, p.n_hill);
    let foxo_on_aid = hill(foxo.clone(), p.half_foxo_on_aid, p.n_hill);

    let pi3k = const_(p.k_pi3k) * tfh.clone();
    let foxo_rhs = const_(p.alpha_f) * help_on_foxo * (const_(1.0) - pi3k.clone())
        - const_(p.delta_f) * foxo.clone()
        + const_(0.08) * bcl6.clone() * foxo.clone() * (const_(1.0) - foxo.clone());

    let bcl6_rhs = const_(p.alpha_b) * foxo.clone() * bcl6.clone() * (const_(1.0) - bcl6.clone())
        + const_(0.12) * foxo.clone() * (const_(1.0) - bcl6.clone())
        - const_(p.delta_b) * bcl6.clone();

    let aid_drive = const_(0.55) * bcl6_on_aid + const_(0.45) * foxo_on_aid;
    let aicda_rhs =
        const_(p.alpha_a) * aid_drive * (const_(1.0) - aicda.clone()) - const_(p.delta_a) * aicda.clone();

    let cxcr4_rhs = const_(p.alpha_x) * cxcl12_on_cxcr4 * foxo.clone() * (const_(1.0) - cxcr4.clone())
        - const_(p.delta_x) * cxcr4.clone();

    vec![
        OdeRule {
            cell_type: "gc_B_cell".into(),
            behavior: "FOXO1".into(),
            rate_expr: foxo_rhs,
            bounds: Some((0.0, 1.0)),
        },
        OdeRule {
            cell_type: "gc_B_cell".into(),
            behavior: "BCL6".into(),
            rate_expr: bcl6_rhs,
            bounds: Some((0.0, 1.0)),
        },
        OdeRule {
            cell_type: "gc_B_cell".into(),
            behavior: "AICDA".into(),
            rate_expr: aicda_rhs,
            bounds: Some((0.0, 1.0)),
        },
        OdeRule {
            cell_type: "gc_B_cell".into(),
            behavior: "CXCR4".into(),
            rate_expr: cxcr4_rhs,
            bounds: Some((0.0, 1.0)),
        },
    ]
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
