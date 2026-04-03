//! Gene matching, synonym lookup, and sequence extraction.

use std::collections::HashMap;

use crate::genbank::{Feature, GenBankRecord};

/// Built-in gene synonym table for mitochondrial genes.
pub fn gene_synonyms() -> HashMap<&'static str, Vec<&'static str>> {
    HashMap::from([
        ("12S", vec!["12S", "12S rRNA", "s-rRNA", "12S ribosomal RNA", "rrnS"]),
        ("16S", vec!["16S", "16S rRNA", "l-rRNA", "16S ribosomal RNA", "rrnL"]),
        ("COI", vec!["COI", "COXI", "COX1", "CO1", "cytochrome oxidase subunit 1", "cytochrome c oxidase subunit I", "cytochrome oxidase subunit I"]),
        ("COII", vec!["COII", "COXII", "COX2", "CO2", "cytochrome oxidase subunit 2", "cytochrome c oxidase subunit II"]),
        ("COIII", vec!["COIII", "COXIII", "COX3", "CO3", "cytochrome oxidase subunit 3"]),
        ("CYTB", vec!["CYTB", "Cyt b", "cytochrome b", "cob", "CYB"]),
        ("ND1", vec!["ND1", "NAD1", "NADH1", "NADH dehydrogenase subunit 1"]),
        ("ND2", vec!["ND2", "NAD2", "NADH2", "NADH dehydrogenase subunit 2"]),
        ("ND3", vec!["ND3", "NAD3", "NADH3", "NADH dehydrogenase subunit 3"]),
        ("ND4", vec!["ND4", "NAD4", "NADH4", "NADH dehydrogenase subunit 4"]),
        ("ND4L", vec!["ND4L", "NAD4L", "NADH4L", "NADH dehydrogenase subunit 4L"]),
        ("ND5", vec!["ND5", "NAD5", "NADH5", "NADH dehydrogenase subunit 5"]),
        ("ND6", vec!["ND6", "NAD6", "NADH6", "NADH dehydrogenase subunit 6"]),
        ("ATP6", vec!["ATP6", "ATPase6", "ATPase 6", "ATP synthase F0 subunit 6"]),
        ("ATP8", vec!["ATP8", "ATPase8", "ATPase 8", "ATP synthase F0 subunit 8"]),
    ])
}

/// Build target gene map from user-specified gene names.
pub fn parse_target_genes(genes: &[String]) -> HashMap<String, Vec<String>> {
    let synonyms = gene_synonyms();
    let mut result = HashMap::new();

    for gene in genes {
        let upper = gene.to_uppercase();
        if let Some(syns) = synonyms.get(upper.as_str()) {
            result.insert(upper, syns.iter().map(|s| s.to_string()).collect());
        } else {
            // Custom gene - basic synonyms
            result.insert(
                upper.clone(),
                vec![gene.clone(), upper.clone(), gene.to_lowercase()],
            );
        }
    }
    result
}

/// Match an annotation name against target genes. Returns the canonical gene name if matched.
pub fn match_gene(annotation: &str, target_genes: &HashMap<String, Vec<String>>) -> Option<String> {
    let ann_lower = annotation.trim().to_lowercase();

    // Exact match first
    for (gene, synonyms) in target_genes {
        for syn in synonyms {
            if ann_lower == syn.to_lowercase() {
                return Some(gene.clone());
            }
        }
    }

    // Partial match (contains)
    for (gene, synonyms) in target_genes {
        for syn in synonyms {
            if ann_lower.contains(&syn.to_lowercase()) {
                return Some(gene.clone());
            }
        }
    }

    None
}

/// Get the gene name from a feature's qualifiers.
pub fn get_gene_name(feature: &Feature) -> Option<&str> {
    feature
        .qualifiers
        .get("gene")
        .or_else(|| feature.qualifiers.get("product"))
        .map(|s| s.as_str())
}

/// Check if a feature is a tRNA.
pub fn is_trna(feature: &Feature, gene_name: &str) -> bool {
    if feature.feature_type == "tRNA" {
        return true;
    }
    let lower = gene_name.to_lowercase();
    lower.starts_with("trn") || lower.contains("trna")
}

/// Extract the DNA sequence for a feature from its parent record.
pub fn extract_sequence(record: &GenBankRecord, feature: &Feature) -> String {
    let seq_bytes = record.sequence.as_bytes();
    let mut result = String::new();

    for loc in &feature.locations {
        let start = loc.start.saturating_sub(1); // 1-indexed to 0-indexed
        let end = loc.end.min(seq_bytes.len());
        if start < seq_bytes.len() {
            result.push_str(&record.sequence[start..end]);
        }
    }

    if feature.complement {
        reverse_complement(&result)
    } else {
        result
    }
}

/// Reverse complement a DNA sequence.
pub fn reverse_complement(seq: &str) -> String {
    seq.bytes()
        .rev()
        .map(|b| match b {
            b'A' => b'T', b'T' => b'A', b'G' => b'C', b'C' => b'G',
            b'a' => b't', b't' => b'a', b'g' => b'c', b'c' => b'g',
            b'N' => b'N', b'n' => b'n',
            b'R' => b'Y', b'Y' => b'R', b'S' => b'S', b'W' => b'W',
            b'K' => b'M', b'M' => b'K', b'B' => b'V', b'V' => b'B',
            b'D' => b'H', b'H' => b'D',
            other => other,
        })
        .map(|b| b as char)
        .collect()
}

/// Clean a name for use in filenames.
pub fn clean_name(name: &str) -> String {
    name.replace(' ', "_")
        .replace('/', "_")
        .replace('\\', "_")
        .replace(':', "_")
}

/// Get a location string for FASTA headers (e.g. "266-1300").
pub fn get_location_string(feature: &Feature) -> String {
    if feature.locations.is_empty() {
        return "1-1".to_string();
    }
    let start = feature.locations.first().unwrap().start;
    let end = feature.locations.last().unwrap().end;
    format!("{}-{}", start, end)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement("ATGC"), "GCAT");
        assert_eq!(reverse_complement("AAAA"), "TTTT");
        assert_eq!(reverse_complement(""), "");
        assert_eq!(reverse_complement("ATNC"), "GNAT");
    }

    #[test]
    fn test_match_gene_exact() {
        let targets = parse_target_genes(&["COI".to_string(), "CYTB".to_string()]);
        assert_eq!(match_gene("COX1", &targets), Some("COI".to_string()));
        assert_eq!(match_gene("cytochrome b", &targets), Some("CYTB".to_string()));
        assert_eq!(match_gene("ND1", &targets), None);
    }

    #[test]
    fn test_match_gene_partial() {
        let targets = parse_target_genes(&["COI".to_string()]);
        assert_eq!(
            match_gene("cytochrome c oxidase subunit I", &targets),
            Some("COI".to_string())
        );
    }

    #[test]
    fn test_is_trna() {
        let feature = Feature {
            feature_type: "tRNA".to_string(),
            locations: vec![],
            complement: false,
            qualifiers: std::collections::HashMap::new(),
        };
        assert!(is_trna(&feature, "tRNA-Leu"));

        let cds = Feature {
            feature_type: "CDS".to_string(),
            locations: vec![],
            complement: false,
            qualifiers: std::collections::HashMap::new(),
        };
        assert!(!is_trna(&cds, "COX1"));
        assert!(is_trna(&cds, "trnL"));
    }

    #[test]
    fn test_clean_name() {
        assert_eq!(clean_name("cytochrome b"), "cytochrome_b");
        assert_eq!(clean_name("tRNA-Leu"), "tRNA-Leu");
    }
}
