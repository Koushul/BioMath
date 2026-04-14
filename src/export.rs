use std::io::{Read, Write};

use csv::WriterBuilder;

use crate::error::BioMathError;
use crate::grammar::rules_from_csv_str;
use crate::model::Model;
use crate::response::ResponseFnKind;
use crate::rule::{Response, Rule};

impl Rule {
    /// Paper-style row: 8 columns (Hill-only) or 11 columns with `response_kind` + two parameters
    /// (see [`crate::grammar`]).
    pub fn to_csv_row(&self) -> String {
        let resp = match self.response {
            Response::Increases => "increases",
            Response::Decreases => "decreases",
        };
        let dead = if self.applies_to_dead { 1 } else { 0 };
        match &self.response_fn {
            ResponseFnKind::Hill {
                half_max,
                hill_power,
            } => format!(
                "{} ; {} ; {} ; {} ; {} ; {} ; {} ; {} ; hill ; {} ; {}",
                self.cell_type,
                self.signal,
                resp,
                self.behavior,
                self.max_response,
                half_max,
                hill_power,
                dead,
                half_max,
                hill_power
            ),
            ResponseFnKind::Linear {
                min_threshold,
                max_threshold,
            } => format!(
                "{} ; {} ; {} ; {} ; {} ; 0 ; 1 ; {} ; linear ; {} ; {}",
                self.cell_type,
                self.signal,
                resp,
                self.behavior,
                self.max_response,
                dead,
                min_threshold,
                max_threshold
            ),
            ResponseFnKind::Step { threshold } => format!(
                "{} ; {} ; {} ; {} ; {} ; 0 ; 1 ; {} ; step ; {} ; 0",
                self.cell_type,
                self.signal,
                resp,
                self.behavior,
                self.max_response,
                dead,
                threshold
            ),
            ResponseFnKind::Custom { .. } => format!(
                "{} ; {} ; {} ; {} ; {} ; 0 ; 1 ; {} ; hill ; 1 ; 1",
                self.cell_type,
                self.signal,
                resp,
                self.behavior,
                self.max_response,
                dead
            ),
        }
    }
}

impl Model {
    pub fn export_summary(&self) -> String {
        let mut s = format!("# Model: {}\n\n", self.name);
        for r in &self.compiled.raw_rules {
            s.push_str(&r.to_english());
            s.push('\n');
        }
        for o in &self.ode_rules {
            s.push_str(&o.to_english());
            s.push('\n');
        }
        s
    }

    pub fn export_html(&self) -> String {
        let mut h = String::from("<!DOCTYPE html><html><head><meta charset=\"utf-8\"><title>");
        h.push_str(&self.name);
        h.push_str("</title></head><body><h1>");
        h.push_str(&self.name);
        h.push_str("</h1><ul>");
        for r in &self.compiled.raw_rules {
            h.push_str("<li><code>");
            h.push_str(&html_escape(&r.to_english()));
            h.push_str("</code></li>");
        }
        for o in &self.ode_rules {
            h.push_str("<li><code>");
            h.push_str(&html_escape(&o.to_latex()));
            h.push_str("</code></li>");
        }
        h.push_str("</ul></body></html>");
        h
    }

    pub fn to_csv_writer(&self, w: impl Write) -> Result<(), BioMathError> {
        let mut wtr = WriterBuilder::new()
            .delimiter(b';')
            .has_headers(false)
            .from_writer(w);
        for r in &self.compiled.raw_rules {
            let fields: Vec<String> = r
                .to_csv_row()
                .split(';')
                .map(|s| s.trim().to_string())
                .collect();
            let fields_ref: Vec<&str> = fields.iter().map(|s| s.as_str()).collect();
            wtr.write_record(fields_ref)
                .map_err(|e| BioMathError::Csv(e.to_string()))?;
        }
        wtr.flush().map_err(|e| BioMathError::Csv(e.to_string()))?;
        Ok(())
    }

    pub fn from_csv_reader(mut r: impl Read) -> Result<Self, BioMathError> {
        let mut buf = String::new();
        r.read_to_string(&mut buf)
            .map_err(|e| BioMathError::Csv(e.to_string()))?;
        let rules = rules_from_csv_str(&buf).map_err(|e| BioMathError::Csv(e))?;
        let compiled = crate::compile::compile(rules, &crate::compile::BaseValueMap::default(), None)
            .map_err(|e| BioMathError::Csv(format!("compile: {e:?}")))?;
        Ok(Model {
            name: "imported".into(),
            compiled,
            ode_rules: Vec::new(),
        })
    }
}

fn html_escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}
