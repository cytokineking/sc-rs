# sc-rs

An open-source Rust implementation of the Shape Complementarity (SC) "goodness of fit" between two protein interfaces as defined by Lawrence & Colman (1993). The SC statistic ranges from 0 to 1, where 1 indicates maximum shape compatibility.

## What is this?
- Rust reimplementation of the classic SC algorithm for protein–protein interfaces
- CPU multi-core by default (Rayon); Connolly-style molecular surface generation (convex, toroidal, concave), peripheral trimming, nearest-neighbor statistics, and SC medians
- Protein-only by default (PDB `ATOM` records); `HETATM` (ions, solvent, ligands) ignored
- CLI for quick runs on PDBs; library API for embedding in other tools

## Quickstart
Build and run the CLI:
```bash
# Debug
cargo run --bin sc -- test-pdb.pdb A B --json

# Release (faster)
cargo run --release --bin sc -- test-pdb.pdb A B --json
```
Example output:
```json
{
  "version": "1.0.0",
  "sc": 0.5602209022718301,
  "median_distance": 0.7034137682650752,
  "trimmed_area": 825.369136270794,
  "atoms_mol1": 1852,
  "atoms_mol2": 1295,
  "elapsed_ms": 1234
}
```

## Installation
You need Rust (1.70+ recommended):
```bash
# macOS (with Homebrew)
brew install rust

# Verify	rustc --version
cargo --version
```

## Project layout
```
src/
  sc/
    mod.rs                         # module wiring & re-exports
    vector3.rs                     # small Vec3 math utility (dot/cross/normalize)
    types.rs                       # Atom, Dot, Probe, Results, enums
    settings.rs                    # defaults and paper-backed constants
    atomic_radii.rs                # JSON radii loader + wildcard matching
    surface_generator.rs           # Connolly surfaces & dot generation
    sc_calculator.rs               # Trimming, NN medians, SC (no histograms)
  lib.rs                           # library entry (exports sc module)

src/bin/sc.rs                      # CLI: PDB parsing, chain split, run SC
atomic_radii.json                  # embedded default atomic radii
```

## CLI usage
```bash
# Show SC in plain text
cargo run --bin sc -- test-pdb.pdb A B

# JSON output (easier to script)
cargo run --bin sc -- test-pdb.pdb A B --json

# Disable parallelization (for benchmarking or debugging)
cargo run --bin sc -- test-pdb.pdb A B --json --no-parallel

# Run the compiled binary directly
# Debug:   target/debug/sc
# Release: target/release/sc
```

## Parallelization
- Parallel processing is enabled by default using Rayon and will automatically use available logical CPUs.
- Disable with the CLI flag `--no-parallel` or in code via `sc.settings_mut().enable_parallel = false;`.
- Control threads with the environment variable `RAYON_NUM_THREADS` (e.g., `RAYON_NUM_THREADS=8`).
- Parallelized stages: peripheral band trimming and nearest-neighbor pairing. Results are deterministic and unaffected by parallelism.

## PDB parsing
- Only standard protein `ATOM` records are loaded; `HETATM` (ions, solvent, ligands) are ignored by default.
- Hydrogens are skipped.
- Future direction: add optional support for additional ligands by extending the atomic radii table with their residue/atom patterns.

## Library usage (embed in your Rust app)
```rust
use sc_rs::sc::types::Atom;
use sc_rs::sc::vector3::Vec3;
use sc_rs::sc::ScCalculator;

fn compute_sc(atoms_a: Vec<(Vec3, String, String)>, atoms_b: Vec<(Vec3, String, String)>) -> anyhow::Result<f64> {
    let mut sc = ScCalculator::new();

    for (pos, atom_name, res_name) in atoms_a {
        let mut a = Atom::new();
        a.coor = pos;
        a.atom = atom_name;
        a.residue = res_name;
        sc.add_atom(0, a)?; // molecule 1
    }
    for (pos, atom_name, res_name) in atoms_b {
        let mut a = Atom::new();
        a.coor = pos;
        a.atom = atom_name;
        a.residue = res_name;
        sc.add_atom(1, a)?; // molecule 2
    }

    let results = sc.calc()?;
    Ok(results.sc)
}
```

## Radii
- Default radii are embedded in the binary at build time from `atomic_radii.json`. You can ship and run the binary without providing any radii file.
- The embedded defaults are selected for maximum compatibility with widely used SC workflows; users may substitute their own radii without rebuilding.
- Optional override: set an environment variable to use a custom JSON file at runtime (no rebuild needed):
  - `ATOMIC_RADII=/path/to/atomic_radii.json`
  - or `ATOMIC_RADII_PATH=/path/to/atomic_radii.json`

### Format (JSON)
- An array of objects: `{ "residue": "GLU", "atom": "OE*", "radius": 1.60 }`
- Generic entries use `***` for the residue to match any residue.
- Atom patterns may end with `*` to mean a prefix match (e.g., `OE*` matches `OE1`, `OE2`).

Examples:
```json
[
  {"residue":"GLU","atom":"OE*","radius":1.60},
  {"residue":"GLU","atom":"CG","radius":1.90},
  {"residue":"***","atom":"N","radius":1.65}
]
```

### Matching semantics
- Exact match when there is no `*` in the pattern.
- With a trailing `*`, only the prefix before `*` must match.
- `***` matches any residue name.
- First match wins: earlier entries take precedence.
- Debug with `ATOMIC_RADII_DEBUG=1`.
- Element fallback: if no explicit pattern matches, the first letter of the atom name is used to try a generic entry (e.g., `***:C`, `***:N`).

Notes:
- To include additional ligands in calculations, add their residue and atom patterns to your radii JSON (or provide a custom file via `ATOMIC_RADII`/`ATOMIC_RADII_PATH`).

## Paper-backed constants and definitions
- Lawrence & Colman (1993), Fig. 1:
  - Gaussian weight w = 0.5 Å^-2
  - Peripheral band d = 1.5 Å excluded from the periphery
  - ~15 dots per Å^2 is sufficient; doubling density does not materially change Sc
- Implementation uses buried surface points per the paper; normals are outward/inward unit normals.
- SA→B = (n_A · n_B) exp(−w |x_A − x_B|^2) with nearest-neighbor mapping of each buried x_A to x_B; Sc is the average of medians of SA→B and SB→A. Medians are used due to skew.

## Connolly surface representation
- Generating probe-based contact and re-entrant surfaces follows Connolly (1983). We use convex, toroidal (re-entrant), and concave patches.
- Stability is achieved with sufficiently high dot density; doubling density does not materially change Sc (Lawrence & Colman, 1993).

## Design notes
- Taken directly from Lawrence & Colman (1993): w=0.5 Å^-2, d=1.5 Å peripheral exclusion, ~15 dots/Å^2 sampling target, medians (not means), NN pairing, Sc formula.
- Implementation choices: sampling method (fixed-step arcs/circles), neighbor acceleration and epsilon, radii source (`atomic_radii.json`), trimming algorithm (k-NN-style erosion), data structures (indices + enums), control flow (collect_neighbors → build_probes → emit_surface_patches), and naming.

## Algorithm notes
- Geometry follows the Connolly-style surface: convex (accessible), toroidal (re-entrant), and concave (probe triangle) patches.
- After peripheral trimming, nearest-neighbor distances and outward normal products are used to compute medians directly (no histograms), with Gaussian weighting exp(−w r^2) using w=0.5 Å^-2.

## Why Rust
- Memory safety and performance for geometry-heavy kernels
- Built-in parallelism via Rayon with deterministic results; auto-uses available CPUs
- Simple distribution: one cross-platform binary; easy CI/CD
- Clean, typed library API plus CLI for scripting and embedding
- Extensible configuration (JSON radii, environment overrides) and reproducible builds

## Contributing
- Issues and PRs welcome. Please run:
```bash
cargo fmt
cargo clippy
cargo test   # (tests to be added)
```
