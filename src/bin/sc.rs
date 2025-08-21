use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};

use sc_rs::sc::types::{Atom, Results};
use sc_rs::sc::vector3::Vec3;
use sc_rs::sc::ScCalculator;

#[derive(serde::Serialize)]
struct Output {
    version: &'static str,
    sc: f64,
    median_distance: f64,
    trimmed_area: f64,
    atoms_mol1: usize,
    atoms_mol2: usize,
    elapsed_ms: u128,
}

fn parse_pdb_atoms(path: &str, chain1: &str, chain2: &str) -> anyhow::Result<(Vec<(Vec3, String, String, String)>, Vec<(Vec3, String, String, String)>)> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut mol1 = Vec::new();
    let mut mol2 = Vec::new();
    for line in reader.lines() {
        let l = line?;
        // Use only standard protein ATOM records; ignore ligands/ions/water in HETATM
        if l.starts_with("ATOM") {
            if l.len() < 54 { continue; }
            // Skip alternate locations other than ' ' or 'A' to mirror common PDB handling
            let alt = if l.len() >= 17 { l[16..17].chars().next().unwrap_or(' ') } else { ' ' };
            if alt != ' ' && alt != 'A' { continue; }
            let atom_name = l[12..16].trim().to_string();
            // Skip hydrogens (use heavy atoms only)
            let element = if l.len() >= 78 { l[76..78].trim().to_string() } else { String::new() };
            if element.eq_ignore_ascii_case("H") || atom_name.starts_with('H') || atom_name.ends_with('H') || atom_name.contains("H") && atom_name.chars().next().unwrap_or(' ').is_ascii_digit() {
                continue;
            }
            let res_name = if l.len() >= 20 { l[17..20].trim().to_string() } else { String::from("UNK") };
            let chain_id = if l.len() >= 22 { l[21..22].to_string() } else { String::from(" ") };
            let x: f64 = l[30..38].trim().parse().unwrap_or(0.0);
            let y: f64 = l[38..46].trim().parse().unwrap_or(0.0);
            let z: f64 = l[46..54].trim().parse().unwrap_or(0.0);
            let rec = (Vec3::new(x,y,z), atom_name, res_name, chain_id.clone());
            if chain_id == chain1 { mol1.push(rec); }
            else if chain_id == chain2 { mol2.push(rec); }
        }
    }
    Ok((mol1, mol2))
}

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = env::args().collect();
    if args.len() < 4 {
        eprintln!("Usage: sc <pdb_file> <chain1> <chain2> [--json] [--no-parallel]");
        std::process::exit(1);
    }
    let pdb = &args[1];
    let chain1 = &args[2];
    let chain2 = &args[3];
    let json = args.iter().any(|a| a == "--json");
    let no_parallel = args.iter().any(|a| a == "--no-parallel");

    let (mol1, mol2) = parse_pdb_atoms(pdb, chain1, chain2)?;
    if mol1.is_empty() || mol2.is_empty() {
        anyhow::bail!("No atoms found for one or both chains");
    }

    let mut sc = ScCalculator::new();
    if no_parallel { sc.settings_mut().enable_parallel = false; }
    // Defaults already set; keep them
    for (pos, atom_name, res_name, _chain) in mol1.iter() {
        let mut a = Atom::new();
        a.coor = *pos;
        a.atom = atom_name.clone();
        a.residue = res_name.clone();
        sc.add_atom(0, a)?;
    }
    for (pos, atom_name, res_name, _chain) in mol2.iter() {
        let mut a = Atom::new();
        a.coor = *pos;
        a.atom = atom_name.clone();
        a.residue = res_name.clone();
        sc.add_atom(1, a)?;
    }

    let t0 = std::time::Instant::now();
    let results: Results = sc.calc()?;
    let elapsed = t0.elapsed().as_millis();
    if json {
        let out = Output { version: env!("CARGO_PKG_VERSION"), sc: results.sc, median_distance: results.distance, trimmed_area: results.area, atoms_mol1: results.surfaces[0].n_atoms, atoms_mol2: results.surfaces[1].n_atoms, elapsed_ms: elapsed };
        println!("{}", serde_json::to_string_pretty(&out)?);
    } else {
        println!("SC: {:.3}", results.sc);
        println!("Median distance: {:.3}", results.distance);
        println!("Trimmed area: {:.3}", results.area);
        println!("Atoms: {} + {}", results.surfaces[0].n_atoms, results.surfaces[1].n_atoms);
        println!("Elapsed: {} ms", elapsed);
    }
    Ok(())
}
