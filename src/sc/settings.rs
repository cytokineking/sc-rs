/// Lawrence & Colman (1993), Fig. 1: Gaussian weight w = 0.5 Å^-2
pub const GAUSSIAN_W: f64 = 0.5;
/// Lawrence & Colman (1993): exclude band d = 1.5 Å from periphery
pub const PERIPH_BAND: f64 = 1.5;
/// Lawrence & Colman (1993): ~15 dots per Å^2 sufficient; doubling density does not materially change Sc
pub const DOT_DENSITY: f64 = 15.0;

#[derive(Clone, Debug)]
pub struct Settings {
	/// Probe radius (Connolly 1983)
	pub rp: f64,
	/// Target dot density per Å^2 (Lawrence & Colman 1993)
	pub dot_density: f64,
	/// Peripheral exclusion band d in Å (Lawrence & Colman 1993)
	pub peripheral_band: f64,
	/// Heuristic separation cutoff for attention classification (implementation choice)
	pub separation_cutoff: f64,
	/// Gaussian weight parameter w in Å^-2 (Lawrence & Colman 1993)
	pub gaussian_w: f64,
	/// Prefer using provided per-atom type radii when available (implementation choice)
	pub use_atom_type_radius: bool,
	/// Enable Rayon-parallel sections (trimming and neighbor pairing)
	pub enable_parallel: bool,
}

impl Default for Settings {
	fn default() -> Self {
		Self {
			rp: 1.7,
			dot_density: DOT_DENSITY,
			peripheral_band: PERIPH_BAND,
			separation_cutoff: 8.0,
			gaussian_w: GAUSSIAN_W,
			use_atom_type_radius: false,
			enable_parallel: true,
		}
	}
}
