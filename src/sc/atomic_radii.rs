use std::fs::File;
use std::io::{self, Read};

use crate::sc::types::{AtomRadius, ScValue};

#[derive(serde::Deserialize)]
struct RadiusRecord { residue: String, atom: String, radius: ScValue }

pub fn read_atomic_radii_from_path(path: &str) -> io::Result<Vec<AtomRadius>> {
	let mut f = File::open(path)?;
	let mut buf = String::new();
	f.read_to_string(&mut buf)?;
	read_atomic_radii_from_str(&buf)
}

pub fn read_atomic_radii_from_str(data: &str) -> io::Result<Vec<AtomRadius>> {
	let recs: Vec<RadiusRecord> = serde_json::from_str(data)
		.map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("invalid radii json: {e}")))?;
	Ok(recs.into_iter().filter(|r| r.radius > 0.0).map(|r| AtomRadius { residue: r.residue, atom: r.atom, radius: r.radius }).collect())
}

pub fn embedded_atomic_radii() -> Vec<AtomRadius> {
	let data: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/atomic_radii.json"));
	read_atomic_radii_from_str(data).unwrap_or_default()
}

pub fn wildcard_match(query: &str, pattern: &str) -> bool {
	fn rtrim_spaces(s: &str) -> &str {
		let mut end = s.len();
		let b = s.as_bytes();
		while end > 0 && (b[end - 1] as char) == ' ' { end -= 1; }
		&s[..end]
	}

	let q = rtrim_spaces(query);
	let p = rtrim_spaces(pattern);

	if p.starts_with('*') { return true; }

	if let Some(star) = p.find('*') {
		let plen = star;
		if q.len() < plen { return false; }
		return q[..plen] == p[..plen];
	}

	// No '*' in pattern: exact match only to avoid unintended generic fallbacks
	q == p
}

