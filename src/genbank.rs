//! GenBank file parser with support for multi-record files and nested locations.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum ParseError {
    #[error("failed to open file: {0}")]
    Io(#[from] std::io::Error),
    #[error("invalid location: {0}")]
    #[allow(dead_code)]
    InvalidLocation(String),
}

#[derive(Debug, Clone)]
pub struct Location {
    pub start: usize,
    pub end: usize,
}

#[derive(Debug, Clone)]
pub struct Feature {
    pub feature_type: String,
    pub locations: Vec<Location>,
    pub complement: bool,
    pub qualifiers: HashMap<String, String>,
}

#[derive(Debug)]
pub struct GenBankRecord {
    pub accession: String,
    pub organism: String,
    pub sequence: String,
    pub features: Vec<Feature>,
}

const VALID_FEATURE_TYPES: &[&str] = &["rRNA", "CDS", "gene", "tRNA", "misc_RNA"];

/// Parse all records from a GenBank file (supports multi-record files).
pub fn parse_genbank(path: &Path) -> Result<Vec<GenBankRecord>, ParseError> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut records = Vec::new();

    let mut accession = String::new();
    let mut organism = String::new();
    let mut features: Vec<Feature> = Vec::new();
    let mut seq = String::new();
    let mut in_features = false;
    let mut in_origin = false;
    let mut current_feature_type = String::new();
    let mut current_feature_lines: Vec<String> = Vec::new();

    for line in reader.lines() {
        let line = line?;

        if line.starts_with("//") {
            // Finalize pending feature
            if !current_feature_type.is_empty() {
                if let Some(f) = parse_feature_block(&current_feature_type, &current_feature_lines)
                {
                    features.push(f);
                }
            }

            // If we have sequence data, emit a record
            if !seq.is_empty() || !features.is_empty() {
                let acc = if accession.is_empty() {
                    path.file_stem()
                        .unwrap_or_default()
                        .to_string_lossy()
                        .into_owned()
                } else {
                    accession.clone()
                };
                records.push(GenBankRecord {
                    accession: acc,
                    organism: organism.clone(),
                    sequence: std::mem::take(&mut seq),
                    features: std::mem::take(&mut features),
                });
            }

            // Reset for next record
            accession.clear();
            organism.clear();
            in_features = false;
            in_origin = false;
            current_feature_type.clear();
            current_feature_lines.clear();
            continue;
        }

        // VERSION line
        if line.starts_with("VERSION") {
            if let Some(ver) = line.split_whitespace().nth(1) {
                accession = ver.to_string();
            }
            continue;
        }

        // ORGANISM
        if line.contains("ORGANISM") && organism.is_empty() {
            if let Some(idx) = line.find("ORGANISM") {
                organism = line[idx + 8..].trim().to_string();
            }
            continue;
        }

        // Section markers
        if line.starts_with("FEATURES") {
            in_features = true;
            continue;
        }
        if line.starts_with("ORIGIN") {
            in_features = false;
            in_origin = true;
            // Save pending feature
            if !current_feature_type.is_empty() {
                if let Some(f) = parse_feature_block(&current_feature_type, &current_feature_lines)
                {
                    features.push(f);
                }
                current_feature_type.clear();
                current_feature_lines.clear();
            }
            continue;
        }

        // Features section
        if in_features {
            // New feature: type starts at column 5, non-space
            if line.len() > 5 && line.as_bytes()[5] != b' ' {
                // Save previous feature
                if !current_feature_type.is_empty() {
                    if let Some(f) =
                        parse_feature_block(&current_feature_type, &current_feature_lines)
                    {
                        features.push(f);
                    }
                }

                let fields: Vec<&str> = line.split_whitespace().collect();
                if fields.len() >= 2 {
                    current_feature_type = fields[0].to_string();
                    current_feature_lines = vec![line.to_string()];
                }
            } else if !current_feature_type.is_empty() {
                current_feature_lines.push(line.to_string());
            }
        }

        // Sequence in ORIGIN
        if in_origin {
            for ch in line.bytes() {
                if ch.is_ascii_alphabetic() {
                    seq.push(ch.to_ascii_uppercase() as char);
                }
            }
        }
    }

    // Handle file without trailing //
    if !seq.is_empty() || !features.is_empty() {
        if !current_feature_type.is_empty() {
            if let Some(f) = parse_feature_block(&current_feature_type, &current_feature_lines) {
                features.push(f);
            }
        }
        let acc = if accession.is_empty() {
            path.file_stem()
                .unwrap_or_default()
                .to_string_lossy()
                .into_owned()
        } else {
            accession
        };
        records.push(GenBankRecord {
            accession: acc,
            organism,
            sequence: seq,
            features,
        });
    }

    Ok(records)
}

/// Parse a feature block into a Feature struct.
fn parse_feature_block(feature_type: &str, lines: &[String]) -> Option<Feature> {
    if !VALID_FEATURE_TYPES.contains(&feature_type) {
        return None;
    }
    if lines.is_empty() {
        return None;
    }

    // First line has the location
    let fields: Vec<&str> = lines[0].split_whitespace().collect();
    if fields.len() < 2 {
        return None;
    }

    let (complement, locations) = parse_location(fields[1]);

    // Parse qualifiers
    let mut qualifiers = HashMap::new();
    let mut current_key = String::new();
    let mut current_val = String::new();

    for line in &lines[1..] {
        let trimmed = line.trim();
        if let Some(rest) = trimmed.strip_prefix('/') {
            // Save previous qualifier
            if !current_key.is_empty() {
                qualifiers.insert(
                    std::mem::take(&mut current_key),
                    current_val.trim_matches('"').to_string(),
                );
                current_val.clear();
            }

            if let Some(eq_idx) = rest.find('=') {
                current_key = rest[..eq_idx].to_string();
                current_val = rest[eq_idx + 1..].trim_start_matches('"').to_string();
            } else {
                current_key = rest.to_string();
            }
        } else if !current_key.is_empty() {
            // Continuation line
            current_val.push_str(trimmed.trim_matches('"'));
        }
    }
    if !current_key.is_empty() {
        qualifiers.insert(current_key, current_val.trim_matches('"').to_string());
    }

    Some(Feature {
        feature_type: feature_type.to_string(),
        locations,
        complement,
        qualifiers,
    })
}

/// Parse a location string, handling complement() and join() nesting.
/// Returns (is_complement, locations).
pub fn parse_location(loc_str: &str) -> (bool, Vec<Location>) {
    let mut s = loc_str;
    let mut complement = false;

    // Strip complement()
    if let Some(inner) = s
        .strip_prefix("complement(")
        .and_then(|r| r.strip_suffix(')'))
    {
        complement = true;
        s = inner;
    }

    // Strip join()
    if let Some(inner) = s.strip_prefix("join(").and_then(|r| r.strip_suffix(')')) {
        s = inner;
    }

    // Also handle complement(join(...)) where complement was already stripped
    // and join() inside complement
    if let Some(inner) = s.strip_prefix("join(").and_then(|r| r.strip_suffix(')')) {
        s = inner;
    }

    let locations = parse_ranges(s);
    (complement, locations)
}

/// Parse comma-separated ranges like "1..100,200..300"
fn parse_ranges(s: &str) -> Vec<Location> {
    let mut locs = Vec::new();
    for part in s.split(',') {
        let part = part.trim();
        // Remove partial indicators < >
        let clean: String = part.chars().filter(|&c| c != '<' && c != '>').collect();

        if let Some((start_s, end_s)) = clean.split_once("..") {
            if let (Ok(start), Ok(end)) = (start_s.parse::<usize>(), end_s.parse::<usize>()) {
                locs.push(Location { start, end });
            }
        } else if let Ok(pos) = clean.parse::<usize>() {
            locs.push(Location {
                start: pos,
                end: pos,
            });
        }
    }
    locs
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_simple_range() {
        let (comp, locs) = parse_location("266..1300");
        assert!(!comp);
        assert_eq!(locs.len(), 1);
        assert_eq!(locs[0].start, 266);
        assert_eq!(locs[0].end, 1300);
    }

    #[test]
    fn test_parse_complement() {
        let (comp, locs) = parse_location("complement(133..201)");
        assert!(comp);
        assert_eq!(locs.len(), 1);
        assert_eq!(locs[0].start, 133);
        assert_eq!(locs[0].end, 201);
    }

    #[test]
    fn test_parse_join() {
        let (comp, locs) = parse_location("join(1..100,200..300)");
        assert!(!comp);
        assert_eq!(locs.len(), 2);
        assert_eq!(locs[0].start, 1);
        assert_eq!(locs[0].end, 100);
        assert_eq!(locs[1].start, 200);
        assert_eq!(locs[1].end, 300);
    }

    #[test]
    fn test_parse_complement_join() {
        let (comp, locs) = parse_location("complement(join(1..100,200..300))");
        assert!(comp);
        assert_eq!(locs.len(), 2);
    }

    #[test]
    fn test_parse_partial_range() {
        let (comp, locs) = parse_location("<1487..3017");
        assert!(!comp);
        assert_eq!(locs.len(), 1);
        assert_eq!(locs[0].start, 1487);
        assert_eq!(locs[0].end, 3017);
    }

    #[test]
    fn test_parse_genbank_file() {
        let path = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .parent()
            .unwrap()
            .join("depricated/go/seqripper/mitogenomes_gb/MK642293.gb");
        if path.exists() {
            let records = parse_genbank(&path).unwrap();
            assert_eq!(records.len(), 1);
            assert_eq!(records[0].accession, "MK642293.1");
            assert_eq!(records[0].organism, "Cinygmina furcata");
            assert!(!records[0].sequence.is_empty());
            assert!(!records[0].features.is_empty());
        }
    }
}
