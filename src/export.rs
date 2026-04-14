use std::io::{Read, Write};

use csv::WriterBuilder;

use crate::error::BioMathError;
use crate::model::Model;
use crate::response::ResponseFnKind;
use crate::rule::{Response, Rule};

impl Rule {
    pub fn to_csv_row(&self) -> String {
        let resp = match self.response {
            Response::Increases => "increases",
            Response::Decreases => "decreases",
        };
        let (half_max, hill_power) = match &self.response_fn {
            ResponseFnKind::Hill {
                half_max,
                hill_power,
            } => (*half_max, *hill_power),
            ResponseFnKind::Linear { .. } | ResponseFnKind::Step { .. } => (0.0, 1.0),
            ResponseFnKind::Custom { .. } => (0.0, 1.0),
        };
        let dead = if self.applies_to_dead { 1 } else { 0 };
        format!(
            "{} ; {} ; {} ; {} ; {} ; {} ; {} ; {}",
            self.cell_type,
            self.signal,
            resp,
            self.behavior,
            self.max_response,
            half_max,
            hill_power,
            dead
        )
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
            let resp = match r.response {
                Response::Increases => "increases",
                Response::Decreases => "decreases",
            };
            let (half_max, hill_power) = match &r.response_fn {
                ResponseFnKind::Hill {
                    half_max,
                    hill_power,
                } => (half_max.to_string(), hill_power.to_string()),
                _ => ("0".into(), "1".into()),
            };
            let max_s = r.max_response.to_string();
            wtr.write_record([
                r.cell_type.as_str(),
                r.signal.as_str(),
                resp,
                r.behavior.as_str(),
                max_s.as_str(),
                half_max.as_str(),
                hill_power.as_str(),
                if r.applies_to_dead { "1" } else { "0" },
            ])
            .map_err(|e| BioMathError::Csv(e.to_string()))?;
        }
        wtr.flush().map_err(|e| BioMathError::Csv(e.to_string()))?;
        Ok(())
    }

    pub fn from_csv_reader(mut r: impl Read) -> Result<Self, BioMathError> {
        let mut buf = String::new();
        r.read_to_string(&mut buf)
            .map_err(|e| BioMathError::Csv(e.to_string()))?;
        let mut rules = Vec::new();
        for line in buf.lines() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let parts: Vec<&str> = line.split(';').map(|p| p.trim()).collect();
            if parts.len() < 8 {
                continue;
            }
            let cell_type = parts[0].to_string();
            let signal = parts[1].to_string();
            let response = match parts[2].to_lowercase().as_str() {
                "increases" => Response::Increases,
                _ => Response::Decreases,
            };
            let behavior = parts[3].to_string();
            let max_response: f64 = parts[4].parse().map_err(|e| {
                BioMathError::Csv(format!("max_response: {e}"))
            })?;
            let half_max: f64 = parts[5].parse().map_err(|e| {
                BioMathError::Csv(format!("half_max: {e}"))
            })?;
            let hill_power: f64 = parts[6].parse().map_err(|e| {
                BioMathError::Csv(format!("hill_power: {e}"))
            })?;
            let applies_to_dead = parts[7] == "1";
            rules.push(Rule {
                cell_type,
                signal,
                response,
                behavior,
                max_response,
                response_fn: ResponseFnKind::Hill {
                    half_max,
                    hill_power,
                },
                applies_to_dead,
                condition: None,
            });
        }
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
