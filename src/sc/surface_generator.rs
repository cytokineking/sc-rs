use std::cmp::Ordering;
use std::f64::consts::PI;
use std::env;

use crate::sc::atomic_radii::{read_atomic_radii_from_path, embedded_atomic_radii, wildcard_match};
use crate::sc::settings::Settings;
use crate::sc::types::*;
use crate::sc::vector3::Vec3;
use rayon::prelude::*;
// Geometry was originally single-threaded; selected phases are parallelized when enabled

#[derive(thiserror::Error, Debug)]
pub enum SurfaceCalculatorError {
	#[error("No atoms defined")] NoAtoms,
	#[error("Index out of bounds")] JumpOutOfBounds,
	#[error("Failed to read radii: {0}")] Io(#[from] std::io::Error),
	#[error("Overlapping atoms detected: {0}")] Coincident(String),
	#[error("Geometric construction invalid (far circle) for atom {0}, neighbor {1}")] ImagFar(i32, i32),
	#[error("Geometric construction invalid (containment) for atom {0}, neighbor {1}")] ImagContain(i32, i32),
	#[error("Invalid local frame for atom {0}, neighbor {1}")] NonPositiveFrame(i32, i32),
	#[error("Sampling limit exceeded")] TooManySubdivisions,
}

pub struct SurfaceGenerator {
	pub settings: Settings,
	radii: Vec<crate::sc::types::AtomRadius>,
	pub(crate) run: RunState,
}

#[derive(Clone, Default)]
pub(crate) struct RunState {
	pub atoms: Vec<Atom>,
	pub probes: Vec<Probe>,
	pub dots: [Vec<Dot>; 2],
	pub trimmed_dots: [Vec<usize>; 2],
	pub results: Results,
	pub radmax: ScValue,
}

impl Default for SurfaceGenerator {
	fn default() -> Self { Self::new() }
}

impl SurfaceGenerator {
	pub fn new() -> Self {
		Self { settings: Settings::default(), radii: Vec::new(), run: RunState::default() }
	}

	pub fn init(&mut self) -> Result<(), SurfaceCalculatorError> {
		if self.radii.is_empty() {
			// Default to embedded radii (portable), allow optional override via env
			self.radii = embedded_atomic_radii();
			if let Ok(path) = env::var("ATOMIC_RADII").or_else(|_| env::var("ATOMIC_RADII_PATH")) {
				if let Ok(r) = read_atomic_radii_from_path(&path) { self.radii = r; }
			}
		}
		Ok(())
	}

	pub fn set_radii(&mut self, radii: Vec<crate::sc::types::AtomRadius>) { self.radii = radii; }

	pub fn reset(&mut self) {
		for a in &mut self.run.atoms { a.neighbor_indices.clear(); a.buried_by_indices.clear(); }
		self.run.atoms.clear();
		self.run.probes.clear();
		self.run.dots[0].clear();
		self.run.dots[1].clear();
		self.run.trimmed_dots[0].clear();
		self.run.trimmed_dots[1].clear();
		self.run.results = Results::default();
	}

	pub fn add_atom(&mut self, molecule: i32, mut atom: Atom) -> Result<(), SurfaceCalculatorError> {
		// Ensure radii are loaded before first assignment
		if self.radii.is_empty() { self.init()?; }
		if atom.radius <= 0.0 { self.assign_atom_radius(&mut atom)?; }
		if atom.radius > 0.0 {
			let mol = if molecule == 1 { 1 } else { 0 } as usize;
			atom.density = self.settings.dot_density;
			atom.molecule = mol;
			atom.natom = (self.run.results.n_atoms + 1) as i32;
			atom.accessible = false;
			self.run.atoms.push(atom);
			self.run.results.surfaces[mol].n_atoms += 1;
			self.run.results.n_atoms += 1;
			Ok(())
		} else {
			Err(SurfaceCalculatorError::Io(std::io::Error::other("Failed to assign atom radius")))
		}
	}

	fn assign_atom_radius(&self, atom: &mut Atom) -> Result<(), SurfaceCalculatorError> {
		if self.settings.use_atom_type_radius {
			if atom.atom_type_radius != 0.0 { atom.radius = atom.atom_type_radius; return Ok(()); }
			return Err(SurfaceCalculatorError::Io(std::io::Error::other("Missing atom_type_radius")));
		}
		let debug = env::var("ATOMIC_RADII_DEBUG").ok().map(|v| v == "1" || v.eq_ignore_ascii_case("true")).unwrap_or(false);
		for radius in &self.radii {
			if !wildcard_match(&atom.residue, &radius.residue) { continue; }
			if !wildcard_match(&atom.atom, &radius.atom) { continue; }
			atom.radius = radius.radius;
			if debug {
				let is_generic = radius.residue.starts_with("***");
				eprintln!(
					"[ATOMIC_RADII_DEBUG] matched {}:{} via pattern {}:{} (generic: {}) => {:.2}",
					atom.residue.trim(), atom.atom.trim(), radius.residue.trim(), radius.atom.trim(), is_generic, radius.radius
				);
			}
			return Ok(());
		}
		// Element fallback: if no specific match, use generic element radii (e.g., ***:C)
		let elem = atom.atom.chars().find(|c| c.is_ascii_alphabetic()).map(|c| c.to_ascii_uppercase()).unwrap_or(' ');
		if elem != ' ' {
			let elem_str = elem.to_string();
			for radius in &self.radii {
				if !radius.residue.trim().starts_with("***") { continue; }
				if radius.atom.trim() != elem_str { continue; }
				atom.radius = radius.radius;
				if debug {
					eprintln!(
						"[ATOMIC_RADII_DEBUG] element fallback {}:{} -> ***:{} => {:.2}",
						atom.residue.trim(), atom.atom.trim(), elem_str, radius.radius
					);
				}
				return Ok(());
			}
		}
		if debug { eprintln!("[ATOMIC_RADII_DEBUG] no match for {}:{}", atom.residue.trim(), atom.atom.trim()); }
		Err(SurfaceCalculatorError::Io(std::io::Error::other(format!("No radius for {}:{}", atom.residue, atom.atom))))
	}

	pub fn assign_attention_numbers(&mut self) {
		// Reset per-surface counters before recomputation
		self.run.results.surfaces[0].n_buried_atoms = 0;
		self.run.results.surfaces[0].n_blocked_atoms = 0;
		self.run.results.surfaces[1].n_buried_atoms = 0;
		self.run.results.surfaces[1].n_blocked_atoms = 0;

		let sep2 = self.settings.separation_cutoff * self.settings.separation_cutoff;
		let atoms_len = self.run.atoms.len();
		// Compute min squared distance to any atom in the other molecule, then set attention
		let snapshot: Vec<(usize, usize, f64)> = (0..atoms_len).map(|i| {
			let a1 = &self.run.atoms[i];
			let mut dist_min2 = f64::INFINITY;
			for j in 0..atoms_len {
				let a2 = &self.run.atoms[j];
				if a2.molecule == a1.molecule { continue; }
				let r2 = a1.distance_squared(a2);
				if r2 < dist_min2 { dist_min2 = r2; }
			}
			(i, a1.molecule, dist_min2)
		}).collect();
		for (i, mol, dist_min2) in snapshot {
			let a1 = &mut self.run.atoms[i];
			if dist_min2 >= sep2 {
				a1.attention = Attention::Far;
				self.run.results.surfaces[mol].n_blocked_atoms += 1;
			} else {
				a1.attention = Attention::Buried;
				self.run.results.surfaces[mol].n_buried_atoms += 1;
			}
		}
	}

	pub fn calc(&mut self) -> Result<(), SurfaceCalculatorError> {
		self.init()?;
		self.run.results.valid = 0;
		if self.run.atoms.is_empty() { return Err(SurfaceCalculatorError::NoAtoms); }
		self.assign_attention_numbers();
		self.generate_molecular_surfaces()?;
		Ok(())
	}

	pub(crate) fn generate_molecular_surfaces(&mut self) -> Result<(), SurfaceCalculatorError> {
		if self.run.atoms.is_empty() { return Err(SurfaceCalculatorError::NoAtoms); }
		self.calc_dots_for_all_atoms()?;
		Ok(())
	}

	fn calc_dots_for_all_atoms(&mut self) -> Result<(), SurfaceCalculatorError> {
		self.run.radmax = 0.0;
		for a in &self.run.atoms { if a.radius > self.run.radmax { self.run.radmax = a.radius; } }
		let atoms_ptrs: Vec<*const Atom> = self.run.atoms.iter().map(|a| a as *const Atom).collect();
		let len = self.run.atoms.len();
		// Phase 1: compute neighbors in parallel to avoid repeated mutable borrows
		if self.settings.enable_parallel { self.compute_neighbors_all_parallel()?; }
		for i in 0..len {
			let att = self.run.atoms[i].attention;
			if matches!(att, Attention::Far) { continue; }
			if !self.settings.enable_parallel { let _ = self.find_neighbors_for_atom_by_index(i, &atoms_ptrs)?; }
			if matches!(self.run.atoms[i].attention, Attention::Far) { continue; }
			if matches!(self.run.atoms[i].attention, Attention::Consider) && self.run.atoms[i].buried_by_indices.is_empty() { continue; }
			self.build_probes(i, &atoms_ptrs)?;
			if !self.settings.enable_parallel && self.run.atoms[i].accessible { self.emit_contact_surface_for_atom(i)?; }
		}
		// Phase 3: contact dot generation in parallel (uses per-atom buffers)
		if self.settings.enable_parallel { self.generate_contact_surface_parallel()?; }
		if self.settings.rp > 0.0 {
			if self.settings.enable_parallel { self.generate_concave_surface_parallel()?; }
			else { self.generate_concave_surface()?; }
		}
		Ok(())
	}

	fn compute_neighbors_all_parallel(&mut self) -> Result<(), SurfaceCalculatorError> {
		let len = self.run.atoms.len();
		let rp = self.settings.rp;
		let radmax = self.run.radmax;
		let _bb2 = (4.0 * radmax + 4.0 * rp).powi(2);
		let atoms: &Vec<Atom> = &self.run.atoms;
		let results: Result<Vec<(Vec<usize>, Vec<usize>, bool)>, SurfaceCalculatorError> = (0..len).into_par_iter().map(|i| {
			let atom1 = &atoms[i];
			let mut neighbor_indices: Vec<usize> = Vec::new();
			let mut buried_by_indices: Vec<usize> = Vec::new();
			// count not used; rely on buried_by_indices length
			for j in 0..len {
				if j == i { continue; }
				let atom2 = &atoms[j];
				if atom1.natom == atom2.natom { continue; }
				let d2 = atom1.distance_squared(atom2);
				if atom1.molecule == atom2.molecule {
					if d2 <= 0.0001 {
						return Err(SurfaceCalculatorError::Coincident(format!(
							"{}:{}:{} @ ({:.3},{:.3},{:.3}) == {}:{}:{} @ ({:.3},{:.3},{:.3})",
							atom1.natom, atom1.residue, atom1.atom, atom1.coor.x, atom1.coor.y, atom1.coor.z,
							atom2.natom, atom2.residue, atom2.atom, atom2.coor.x, atom2.coor.y, atom2.coor.z
						)));
					}
					let bridge = atom1.radius + atom2.radius + 2.0 * rp;
					if d2 < bridge * bridge { neighbor_indices.push(j); }
				} else {
					let bridge = atom1.radius + atom2.radius + 2.0 * rp;
					if d2 < bridge * bridge { buried_by_indices.push(j); }
				}
			}
			let center = atom1.coor;
			neighbor_indices.sort_unstable_by(|&a1, &a2| {
				let d1 = atoms[a1].coor.distance_squared(center);
				let d2 = atoms[a2].coor.distance_squared(center);
				if d1 < d2 { Ordering::Less } else if d1 > d2 { Ordering::Greater } else { Ordering::Equal }
			});
			let accessible = neighbor_indices.is_empty();
			Ok((neighbor_indices, buried_by_indices, accessible))
		}).collect();
		let outs = results?;
		for (i, (neighbors, buried_by, accessible)) in outs.into_iter().enumerate() {
			let a1 = &mut self.run.atoms[i];
			a1.neighbor_indices = neighbors;
			a1.buried_by_indices = buried_by;
			if accessible { a1.accessible = true; }
		}
		Ok(())
	}

	fn generate_contact_surface_parallel(&mut self) -> Result<(), SurfaceCalculatorError> {
		let rp = self.settings.rp;
		let atoms: &Vec<Atom> = &self.run.atoms;
		let results: Vec<(usize, Vec<Dot>, usize)> = (0..atoms.len()).into_par_iter().filter_map(|i| {
			let a_i = &atoms[i];
			let att = a_i.attention;
			if matches!(att, Attention::Far) { return None; }
			if matches!(att, Attention::Consider) && a_i.buried_by_indices.is_empty() { return None; }
			if !a_i.accessible { return None; }
			let neighbors = a_i.neighbor_indices.clone();
			let mut north_dir = Vec3::new(0.0, 0.0, 1.0);
			let mut south_dir = Vec3::new(0.0, 0.0, -1.0);
			let mut equatorial_vector = Vec3::new(1.0, 0.0, 0.0);
			let radius_i = a_i.radius;
			let expanded_radius_i = a_i.radius + rp;
			if !neighbors.is_empty() {
				let neighbor = &atoms[neighbors[0]];
				north_dir = a_i.coor - neighbor.coor;
				north_dir.normalize();
				let mut temp_vec = Vec3::new(north_dir.y*north_dir.y + north_dir.z*north_dir.z, north_dir.x*north_dir.x + north_dir.z*north_dir.z, north_dir.x*north_dir.x + north_dir.y*north_dir.y);
				temp_vec.normalize();
				let dt = temp_vec.dot(north_dir);
				if dt.abs() > 0.99 { temp_vec = Vec3::new(1.0, 0.0, 0.0); }
				equatorial_vector = north_dir.cross(temp_vec);
				equatorial_vector.normalize();
				let _ = equatorial_vector.cross(north_dir);
				let radius_neighbor = neighbor.radius;
				let expanded_radius_j = neighbor.radius + rp;
				let dij = a_i.coor.distance(neighbor.coor);
				let unit_axis = (neighbor.coor - a_i.coor) / dij;
				let asymmetry_term = (expanded_radius_i*expanded_radius_i - expanded_radius_j*expanded_radius_j) / dij;
				let midplane_center = (a_i.coor + neighbor.coor) * 0.5 + (unit_axis * (asymmetry_term*0.5));
				let mut far_term = (expanded_radius_i + expanded_radius_j)*(expanded_radius_i + expanded_radius_j) - dij*dij;
				if far_term <= 0.0 { return None; }
				far_term = far_term.sqrt();
				let mut contain_term = dij*dij - (radius_i - radius_neighbor).powi(2);
				if contain_term <= 0.0 { return None; }
				contain_term = contain_term.sqrt();
				let ring_radius = 0.5 * far_term * contain_term / dij;
				let ring_point = midplane_center + (equatorial_vector.cross(north_dir) * ring_radius);
				south_dir = (ring_point - a_i.coor) / expanded_radius_i;
				if north_dir.cross(south_dir).dot(equatorial_vector) <= 0.0 { return None; }
			}
			let mut lats: Vec<Vec3> = Vec::new();
			let o = Vec3::zero();
			let cs = geom_sample_arc(o, radius_i, equatorial_vector, a_i.density, north_dir, south_dir, &mut lats).ok()?;
			if lats.is_empty() { return None; }
			let mut dots: Vec<Dot> = Vec::new();
			let mut points: Vec<Vec3> = Vec::new();
			for ilat in lats.iter() {
				let dt = ilat.dot(north_dir);
				let cen = a_i.coor + (north_dir * dt);
				let mut rad = radius_i*radius_i - dt*dt;
				if rad <= 0.0 { continue; }
				rad = rad.sqrt();
				points.clear();
				let ps = geom_sample_circle(cen, rad, north_dir, a_i.density, &mut points).ok()?;
				if points.is_empty() { continue; }
				let area = ps * cs;
				for &point in points.iter() {
					let pcen = a_i.coor + ((point - a_i.coor) * (expanded_radius_i/radius_i));
					// collision with same-molecule neighbors (skip first neighbor)
					let mut coll = false;
					for &idx in neighbors.iter().skip(1) {
						let a = &atoms[idx];
						if pcen.distance(a.coor) <= (a.radius + rp) { coll = true; break; }
					}
					if coll { continue; }
					// burial check against opposite molecule
					let other_mol = if a_i.molecule == 0 { 1 } else { 0 };
					let mut buried = false;
					for b in atoms.iter() {
						if b.molecule != other_mol { continue; }
						let erl = b.radius + rp;
						let d = pcen.distance_squared(b.coor);
						if d <= erl*erl { buried = true; break; }
					}
					let outnml = if rp <= 0.0 { point - a_i.coor } else { (pcen - point) / rp };
					dots.push(Dot { coor: point, outnml, area, buried, kind: DotKind::Contact, atom_index: i });
				}
			}
			if dots.is_empty() { None } else { let n = dots.len(); Some((a_i.molecule, dots, n)) }
		}).collect();
		for (mol, mut dots, n) in results.into_iter() {
			self.run.results.dots.convex += n;
			self.run.dots[mol].append(&mut dots);
		}
		Ok(())
	}


	fn find_neighbors_for_atom_by_index(&mut self, atom_index: usize, atoms_ptrs: &[*const Atom]) -> Result<bool, SurfaceCalculatorError> {
		let mut nbb = 0;
		let bb2 = (4.0 * self.run.radmax + 4.0 * self.settings.rp).powi(2);
		let total = self.run.atoms.len();
		let (_left, rest) = self.run.atoms.split_at_mut(atom_index);
		let (atom1, _right) = rest.split_first_mut().unwrap();
		atom1.neighbor_indices.clear();
		atom1.buried_by_indices.clear();
		for j in 0..total {
			if j == atom_index { continue; }
			let ptr2 = atoms_ptrs[j];
			let atom2 = unsafe { &*ptr2 };
			if atom1.natom == atom2.natom { continue; }
			if atom1.molecule == atom2.molecule {
				let d2 = atom1.distance_squared(atom2);
				if d2 <= 0.0001 {
					return Err(SurfaceCalculatorError::Coincident(format!(
						"{}:{}:{} @ ({:.3},{:.3},{:.3}) == {}:{}:{} @ ({:.3},{:.3},{:.3})",
						atom1.natom, atom1.residue, atom1.atom, atom1.coor.x, atom1.coor.y, atom1.coor.z,
						atom2.natom, atom2.residue, atom2.atom, atom2.coor.x, atom2.coor.y, atom2.coor.z
					)));
				}
				let bridge = atom1.radius + atom2.radius + 2.0 * self.settings.rp;
				if d2 < bridge * bridge { atom1.neighbor_indices.push(j); }
			} else {
				// Include all opposite-molecule atoms for burial check; geometry will decide actual burial
				let d2 = atom1.distance_squared(atom2);
				if d2 < bb2 { nbb += 1; }
				let bridge = atom1.radius + atom2.radius + 2.0 * self.settings.rp;
				if d2 < bridge * bridge { atom1.buried_by_indices.push(j); }
			}
		}
		if matches!(atom1.attention, Attention::Consider) && nbb == 0 { return Ok(false); }
		if atom1.neighbor_indices.is_empty() { atom1.accessible = true; return Ok(false); }
		let center = atom1.coor;
		atom1.neighbor_indices.sort_unstable_by(|&a1, &a2| {
			let d1 = unsafe { (*atoms_ptrs[a1]).coor.distance_squared(center) };
			let d2 = unsafe { (*atoms_ptrs[a2]).coor.distance_squared(center) };
			if d1 < d2 { Ordering::Less } else if d1 > d2 { Ordering::Greater } else { Ordering::Equal }
		});
		Ok(true)
	}

	fn build_probes(&mut self, atom_index: usize, atoms_ptrs: &[*const Atom]) -> Result<(), SurfaceCalculatorError> {
		let expanded_radius_i;
		let neighbor_indices: Vec<usize>;
		{
			let atom1 = &self.run.atoms[atom_index];
			expanded_radius_i = atom1.radius + self.settings.rp;
			neighbor_indices = atom1.neighbor_indices.clone();
		}
		for &j in &neighbor_indices {
			let atom2 = unsafe { &*atoms_ptrs[j] };
			if atom2.natom <= self.run.atoms[atom_index].natom { continue; }
			let expanded_radius_j = atom2.radius + self.settings.rp;
			let dist_ij = self.run.atoms[atom_index].coor.distance(atom2.coor);
			let unit_axis = (atom2.coor - self.run.atoms[atom_index].coor) / dist_ij;
			let asymmetry_term = (expanded_radius_i*expanded_radius_i - expanded_radius_j*expanded_radius_j) / dist_ij;
			let midplane_center = (self.run.atoms[atom_index].coor + atom2.coor) * 0.5 + (unit_axis * (asymmetry_term*0.5));
			let mut far_term = (expanded_radius_i + expanded_radius_j)*(expanded_radius_i + expanded_radius_j) - dist_ij*dist_ij;
			if far_term <= 0.0 { continue; }
			far_term = far_term.sqrt();
			let mut contain_term = dist_ij*dist_ij - (self.run.atoms[atom_index].radius - atom2.radius).powi(2);
			if contain_term <= 0.0 { continue; }
			contain_term = contain_term.sqrt();
			let ring_radius = 0.5 * far_term * contain_term / dist_ij;
			if neighbor_indices.len() <= 1 {
				self.run.atoms[atom_index].accessible = true;
				self.run.atoms[j].accessible = true;
				break;
			}
			self.build_probe_triplets(atom_index, atoms_ptrs[j], unit_axis, midplane_center, ring_radius)?;
			let has_point_cusp = asymmetry_term.abs() < dist_ij;
			if !matches!(self.run.atoms[atom_index].attention, Attention::Far) || (!matches!(atom2.attention, Attention::Far) && self.settings.rp > 0.0) {
				self.emit_reentrant_surface(atom_index, atoms_ptrs[j], unit_axis, midplane_center, ring_radius, has_point_cusp)?;
			}
		}
		Ok(())
	}

	fn build_probe_triplets(&mut self, atom1_index: usize, atom2_ptr: *const Atom, unit_axis: Vec3, midplane_center: Vec3, ring_radius: ScValue) -> Result<(), SurfaceCalculatorError> {
		let neighbor_indices = self.run.atoms[atom1_index].neighbor_indices.clone();
		let expanded_radius_i = self.run.atoms[atom1_index].radius + self.settings.rp;
		let atom2 = unsafe { &*atom2_ptr };
		let expanded_radius_j = atom2.radius + self.settings.rp;
		let mut made_probe = false;
		for &k in &neighbor_indices {
			let atom3 = &self.run.atoms[k];
			if atom3.natom <= atom2.natom { continue; }
			let expanded_radius_k = atom3.radius + self.settings.rp;
			let dist_jk = atom2.coor.distance(atom3.coor);
			if dist_jk >= expanded_radius_j + expanded_radius_k { continue; }
			let dist_ik = self.run.atoms[atom1_index].coor.distance(atom3.coor);
			if dist_ik >= expanded_radius_i + expanded_radius_k { continue; }
			if matches!(self.run.atoms[atom1_index].attention, Attention::Far) && matches!(atom2.attention, Attention::Far) && matches!(atom3.attention, Attention::Far) { continue; }
			let unit_axis_ik = (atom3.coor - self.run.atoms[atom1_index].coor) / dist_ik;
			let wedge_angle = unit_axis.dot(unit_axis_ik).acos();
			let sin_wedge = wedge_angle.sin();
			if sin_wedge <= 0.0 { let dtijk2 = midplane_center.distance(atom3.coor); let rkp2 = expanded_radius_k*expanded_radius_k - ring_radius*ring_radius; if dtijk2 < rkp2 { return Ok(()); } continue; }
			let axis_normal = unit_axis.cross(unit_axis_ik) / sin_wedge;
			let perp_tangent = axis_normal.cross(unit_axis);
			let asymmetry_term_ik = (expanded_radius_i*expanded_radius_i - expanded_radius_k*expanded_radius_k) / dist_ik;
			let midpoint_ik = (self.run.atoms[atom1_index].coor + atom3.coor)*0.5 + unit_axis_ik * (asymmetry_term_ik*0.5);
			let mut componentwise = midpoint_ik - midplane_center;
			componentwise = Vec3::new(unit_axis_ik.x * componentwise.x, unit_axis_ik.y * componentwise.y, unit_axis_ik.z * componentwise.z);
			let component_sum = componentwise.x + componentwise.y + componentwise.z;
			let torus_center = midplane_center + perp_tangent * (component_sum / sin_wedge);
			let mut height = expanded_radius_i*expanded_radius_i - torus_center.distance_squared(self.run.atoms[atom1_index].coor);
			if height <= 0.0 { continue; }
			height = height.sqrt();
			for is0 in 1..=2 {
				let sign_choice = 3 - 2*is0;
				let probe_center = torus_center + axis_normal * (height * (sign_choice as f64));
				if self.check_atom_collision2_idx(probe_center, atom2, atom3, &neighbor_indices) { continue; }
				let mut probe = Probe { atom_indices: [0; 3], height, point: probe_center, alt: axis_normal * (sign_choice as f64) };
				if sign_choice > 0 { probe.atom_indices = [atom1_index, atom2.natom as usize - 1, k]; }
				else { probe.atom_indices = [atom2.natom as usize - 1, atom1_index, k]; }
				self.run.probes.push(probe);
				made_probe = true;
			}
		}
		if made_probe { self.run.atoms[atom1_index].accessible = true; }
		Ok(())
	}

	fn emit_reentrant_surface(&mut self, atom1_index: usize, atom2_ptr: *const Atom, unit_axis: Vec3, midplane_center: Vec3, ring_radius: ScValue, has_point_cusp: bool) -> Result<(), SurfaceCalculatorError> {
		let neighbors = self.run.atoms[atom1_index].neighbor_indices.clone();
		let density = (self.run.atoms[atom1_index].density + unsafe { &*atom2_ptr }.density) / 2.0;
		let expanded_radius_i = self.run.atoms[atom1_index].radius + self.settings.rp;
		let expanded_radius_j = unsafe { &*atom2_ptr }.radius + self.settings.rp;
		let roll_circle_radius_i = ring_radius * self.run.atoms[atom1_index].radius / expanded_radius_i;
		let roll_circle_radius_j = ring_radius * unsafe { &*atom2_ptr }.radius / expanded_radius_j;
		let mut belt_radius = ring_radius - self.settings.rp; if belt_radius <= 0.0 { belt_radius = 0.0; }
		let mean_radius = (roll_circle_radius_i + 2.0*belt_radius + roll_circle_radius_j) / 4.0;
		let eccentricity = mean_radius / ring_radius;
		let effective_density = eccentricity*eccentricity*density;
		let mut subs: Vec<Vec3> = Vec::new();
		let ts = self.sample_circle(midplane_center, ring_radius, unit_axis, effective_density, &mut subs)?;
		if subs.is_empty() { return Ok(()) }
		for sub in subs.into_iter() {
			let mut tooclose = false;
			for &ni in &neighbors {
				let neighbor = &self.run.atoms[ni];
				if neighbor.natom == unsafe { &*atom2_ptr }.natom { continue; }
				let expanded_neighbor_radius = neighbor.radius + self.settings.rp;
				let d2 = sub.distance_squared(neighbor.coor);
				if d2 < expanded_neighbor_radius*expanded_neighbor_radius { tooclose = true; break; }
			}
			if tooclose { continue; }
			let ring_point = sub;
			self.run.atoms[atom1_index].accessible = true;
			unsafe { (*(atom2_ptr as *mut Atom)).accessible = true; }
			let vec_pi = (self.run.atoms[atom1_index].coor - ring_point) / expanded_radius_i;
			let vec_pj = (unsafe { &*atom2_ptr }.coor - ring_point) / expanded_radius_j;
			let mut toroid_axis = vec_pi.cross(vec_pj); toroid_axis.normalize();
			let mut cusp_term = self.settings.rp*self.settings.rp - ring_radius*ring_radius;
			let has_cusp_point = cusp_term > 0.0 && has_point_cusp;
			let (arc_end_i, arc_end_j) = if has_cusp_point {
				cusp_term = cusp_term.sqrt();
				let qij = midplane_center - unit_axis * cusp_term;
				let _qjk = midplane_center + unit_axis * cusp_term;
				(((qij - ring_point)/self.settings.rp), Vec3::zero())
			} else {
				let mut pq = vec_pi + vec_pj; pq.normalize();
				(pq, pq)
			};
			let mut dot_tmp = arc_end_i.dot(vec_pi);
			if dot_tmp >= 1.0 || dot_tmp <= -1.0 { return Ok(()); }
			dot_tmp = arc_end_j.dot(vec_pj);
			if dot_tmp >= 1.0 || dot_tmp <= -1.0 { return Ok(()) ; }
			if !matches!(self.run.atoms[atom1_index].attention, Attention::Far) {
				let mut points: Vec<Vec3> = Vec::new();
				let ps = self.sample_arc(ring_point, self.settings.rp, toroid_axis, density, vec_pi, arc_end_i, &mut points)?;
				for &point in points.iter() { let area = ps * ts * self.distance_point_to_line(midplane_center, unit_axis, point) / ring_radius; self.run.results.dots.toroidal += 1; let molecule = self.run.atoms[atom1_index].molecule; self.add_dot(molecule, DotKind::Reentrant, point, area, ring_point, atom1_index); }
			}
			let atom2_attention = unsafe { (*atom2_ptr).attention };
			if !matches!(atom2_attention, Attention::Far) {
				let mut points: Vec<Vec3> = Vec::new();
				let ps = self.sample_arc(ring_point, self.settings.rp, toroid_axis, density, arc_end_j, vec_pj, &mut points)?;
				let atom2_index = unsafe { &*atom2_ptr }.natom as usize - 1;
				for &point in points.iter() { let area = ps * ts * self.distance_point_to_line(midplane_center, unit_axis, point) / ring_radius; self.run.results.dots.toroidal += 1; let molecule2 = self.run.atoms[atom2_index].molecule; self.add_dot(molecule2, DotKind::Reentrant, point, area, ring_point, atom2_index); }
			}
		}
		Ok(())
	}

	fn emit_contact_surface_for_atom(&mut self, atom_index: usize) -> Result<(), SurfaceCalculatorError> {
		let neighbors = self.run.atoms[atom_index].neighbor_indices.clone();
		let mut north_dir = Vec3::new(0.0, 0.0, 1.0);
		let mut south_dir = Vec3::new(0.0, 0.0, -1.0);
		let mut equatorial_vector = Vec3::new(1.0, 0.0, 0.0);
		let radius_i = self.run.atoms[atom_index].radius;
		let expanded_radius_i = self.run.atoms[atom_index].radius + self.settings.rp;
		if !neighbors.is_empty() {
			let neighbor = &self.run.atoms[neighbors[0]];
			north_dir = self.run.atoms[atom_index].coor - neighbor.coor;
			north_dir.normalize();
			let mut temp_vec = Vec3::new(north_dir.y*north_dir.y + north_dir.z*north_dir.z, north_dir.x*north_dir.x + north_dir.z*north_dir.z, north_dir.x*north_dir.x + north_dir.y*north_dir.y);
			temp_vec.normalize();
			let dt = temp_vec.dot(north_dir);
			if dt.abs() > 0.99 { temp_vec = Vec3::new(1.0, 0.0, 0.0); }
			equatorial_vector = north_dir.cross(temp_vec);
			equatorial_vector.normalize();
			let _ = equatorial_vector.cross(north_dir);
			let radius_neighbor = neighbor.radius;
			let expanded_radius_j = neighbor.radius + self.settings.rp;
			let dij = self.run.atoms[atom_index].coor.distance(neighbor.coor);
			let unit_axis = (neighbor.coor - self.run.atoms[atom_index].coor) / dij;
			let asymmetry_term = (expanded_radius_i*expanded_radius_i - expanded_radius_j*expanded_radius_j) / dij;
			let midplane_center = (self.run.atoms[atom_index].coor + neighbor.coor) * 0.5 + (unit_axis * (asymmetry_term*0.5));
			let mut far_term = (expanded_radius_i + expanded_radius_j)*(expanded_radius_i + expanded_radius_j) - dij*dij;
			if far_term <= 0.0 { return Err(SurfaceCalculatorError::ImagFar(self.run.atoms[atom_index].natom, neighbor.natom)); }
			far_term = far_term.sqrt();
			let mut contain_term = dij*dij - (radius_i - radius_neighbor).powi(2);
			if contain_term <= 0.0 { return Err(SurfaceCalculatorError::ImagContain(self.run.atoms[atom_index].natom, neighbor.natom)); }
			contain_term = contain_term.sqrt();
			let ring_radius = 0.5 * far_term * contain_term / dij;
			let ring_point = midplane_center + (equatorial_vector.cross(north_dir) * ring_radius);
			south_dir = (ring_point - self.run.atoms[atom_index].coor) / expanded_radius_i;
			if north_dir.cross(south_dir).dot(equatorial_vector) <= 0.0 { return Err(SurfaceCalculatorError::NonPositiveFrame(self.run.atoms[atom_index].natom, neighbor.natom)); }
		}
		let mut lats: Vec<Vec3> = Vec::new();
		let o = Vec3::zero();
		let cs = self.sample_arc(o, radius_i, equatorial_vector, self.run.atoms[atom_index].density, north_dir, south_dir, &mut lats)?;
		if lats.is_empty() { return Ok(()); }
		let mut points: Vec<Vec3> = Vec::new();
		for ilat in lats.iter() {
			let dt = ilat.dot(north_dir);
			let cen = self.run.atoms[atom_index].coor + (north_dir * dt);
			let mut rad = radius_i*radius_i - dt*dt;
			if rad <= 0.0 { continue; }
			rad = rad.sqrt();
			points.clear();
			let ps = self.sample_circle(cen, rad, north_dir, self.run.atoms[atom_index].density, &mut points)?;
			if points.is_empty() { continue; }
			let area = ps * cs;
			for &point in points.iter() {
				let pcen = self.run.atoms[atom_index].coor + ((point - self.run.atoms[atom_index].coor) * (expanded_radius_i/radius_i));
				if self.check_point_collision(pcen, &neighbors) { continue; }
				self.run.results.dots.convex += 1;
				let molecule = self.run.atoms[atom_index].molecule;
				self.add_dot(molecule, DotKind::Contact, point, area, pcen, atom_index);
			}
		}
		Ok(())
	}

	fn check_atom_collision2_idx(&self, probe_center: Vec3, atom1: &Atom, atom2: &Atom, neighbor_indices: &Vec<usize>) -> bool {
		for &ni in neighbor_indices {
			let neighbor = &self.run.atoms[ni];
			if neighbor.natom == atom1.natom || neighbor.natom == atom2.natom { continue; }
			if probe_center.distance_squared(neighbor.coor) <= (neighbor.radius + self.settings.rp).powi(2) { return true; }
		}
		false
	}

	fn generate_concave_surface(&mut self) -> Result<(), SurfaceCalculatorError> {
		let mut lowprobs: Vec<usize> = Vec::new();
		for (idx, probe) in self.run.probes.iter().enumerate() { if probe.height < self.settings.rp { lowprobs.push(idx); } }
		for i in 0..self.run.probes.len() {
			let probe = &self.run.probes[i];
			let aidx = probe.atom_indices;
			if matches!(self.run.atoms[aidx[0]].attention, Attention::Consider) && matches!(self.run.atoms[aidx[1]].attention, Attention::Consider) && matches!(self.run.atoms[aidx[2]].attention, Attention::Consider) { continue; }
			let pijk = probe.point; let uijk = probe.alt; let hijk = probe.height; let density = (self.run.atoms[aidx[0]].density + self.run.atoms[aidx[1]].density + self.run.atoms[aidx[2]].density) / 3.0;
			let mut nears: Vec<usize> = Vec::new();
			for &lp in &lowprobs { if lp == i { continue; } let d2 = pijk.distance_squared(self.run.probes[lp].point); if d2 <= 4.0 * self.settings.rp*self.settings.rp { nears.push(lp); } }
			let mut vp = [Vec3::zero();3];
			for k in 0..3 { vp[k] = self.run.atoms[aidx[k]].coor - pijk; vp[k].normalize(); }
			let mut vectors = [Vec3::zero();3];
			vectors[0] = vp[0].cross(vp[1]).normalized();
			vectors[1] = vp[1].cross(vp[2]).normalized();
			vectors[2] = vp[2].cross(vp[0]).normalized();
			let mut dm = -1.0; let mut mm = 0usize;
			for k in 0..3 { let dt = uijk.dot(vp[k]); if dt > dm { dm = dt; mm = k; } }
			let south_dir = uijk * -1.0; let mut arc_axis = vp[mm].cross(south_dir); arc_axis.normalize();
			let mut lats: Vec<Vec3> = Vec::new(); let o = Vec3::zero();
			let cs = self.sample_arc(o, self.settings.rp, arc_axis, density, vp[mm], south_dir, &mut lats)?; if lats.is_empty() { continue; }
			let mut points: Vec<Vec3> = Vec::new();
			for ilat in lats.iter() {
				let dt = ilat.dot(south_dir); let cen = south_dir * dt; let mut rad = self.settings.rp*self.settings.rp - dt*dt; if rad <= 0.0 { continue; } rad = rad.sqrt();
				points.clear(); let ps = self.sample_circle(cen, rad, south_dir, density, &mut points)?; if points.is_empty() { continue; }
				let area = ps * cs;
				for &point in points.iter() {
					let mut bail = false; for v in vectors.iter() { let dt2 = point.dot(*v); if dt2 >= 0.0 { bail = true; break; } } if bail { continue; }
					let point = point + pijk;
					if (hijk < self.settings.rp && !nears.is_empty()) && self.check_probe_collision_idx(point, &nears, self.settings.rp*self.settings.rp) { continue; }
					let mut mc = 0usize; let mut dmin = 2.0 * self.settings.rp; for kk in 0..3 { let d = point.distance(self.run.atoms[aidx[kk]].coor) - self.run.atoms[aidx[kk]].radius; if d < dmin { dmin = d; mc = kk; } }
					let atom_index = aidx[mc]; let molecule = self.run.atoms[atom_index].molecule; self.run.results.dots.concave += 1; self.add_dot(molecule, DotKind::Cavity, point, area, pijk, atom_index);
				}
			}
		}
		Ok(())
	}

	fn generate_concave_surface_parallel(&mut self) -> Result<(), SurfaceCalculatorError> {
		let rp = self.settings.rp;
		let rp2 = rp*rp;
		let atoms: &Vec<Atom> = &self.run.atoms;
		let probes: &Vec<Probe> = &self.run.probes;
		if probes.is_empty() { return Ok(()); }
		let mut lowprobs: Vec<usize> = Vec::new();
		for (idx, probe) in probes.iter().enumerate() { if probe.height < rp { lowprobs.push(idx); } }
		let results: Vec<(Vec<Dot>, Vec<Dot>, usize)> = (0..probes.len()).into_par_iter().filter_map(|i| {
			let probe = &probes[i];
			let aidx = probe.atom_indices;
			// skip if all 3 atoms are Consider
			if matches!(atoms[aidx[0]].attention, Attention::Consider) && matches!(atoms[aidx[1]].attention, Attention::Consider) && matches!(atoms[aidx[2]].attention, Attention::Consider) { return None; }
			let pijk = probe.point; let uijk = probe.alt; let hijk = probe.height;
			let density = (atoms[aidx[0]].density + atoms[aidx[1]].density + atoms[aidx[2]].density) / 3.0;
			// build nears
			let mut nears: Vec<usize> = Vec::new();
			for &lp in &lowprobs { if lp == i { continue; } let d2 = pijk.distance_squared(probes[lp].point); if d2 <= 4.0 * rp2 { nears.push(lp); } }
			let mut vp = [Vec3::zero();3];
			for k in 0..3 { vp[k] = atoms[aidx[k]].coor - pijk; vp[k].normalize(); }
			let mut vectors = [Vec3::zero();3];
			vectors[0] = vp[0].cross(vp[1]).normalized();
			vectors[1] = vp[1].cross(vp[2]).normalized();
			vectors[2] = vp[2].cross(vp[0]).normalized();
			let mut dm = -1.0; let mut mm = 0usize;
			for k in 0..3 { let dt = uijk.dot(vp[k]); if dt > dm { dm = dt; mm = k; } }
			let south_dir = uijk * -1.0; let mut arc_axis = vp[mm].cross(south_dir); arc_axis.normalize();
			let mut lats: Vec<Vec3> = Vec::new(); let o = Vec3::zero();
			let cs = geom_sample_arc(o, rp, arc_axis, density, vp[mm], south_dir, &mut lats).ok()?; if lats.is_empty() { return None; }
			let mut d0: Vec<Dot> = Vec::new();
			let mut d1: Vec<Dot> = Vec::new();
			let mut points: Vec<Vec3> = Vec::new();
			for ilat in lats.iter() {
				let dt = ilat.dot(south_dir); let cen = south_dir * dt; let mut rad = rp2 - dt*dt; if rad <= 0.0 { continue; } rad = rad.sqrt();
				points.clear(); let ps = geom_sample_circle(cen, rad, south_dir, density, &mut points).ok()?; if points.is_empty() { continue; }
				let area = ps * cs;
				for &point in points.iter() {
					let mut bail = false; for v in vectors.iter() { let dt2 = point.dot(*v); if dt2 >= 0.0 { bail = true; break; } } if bail { continue; }
					let point = point + pijk;
					if hijk < rp && !nears.is_empty() {
						let mut coll = false; for &np in &nears { let p = &probes[np]; if point.distance_squared(p.point) < rp2 { coll = true; break; } }
						if coll { continue; }
					}
					let mut mc = 0usize; let mut dmin = 2.0 * rp; for kk in 0..3 { let d = point.distance(atoms[aidx[kk]].coor) - atoms[aidx[kk]].radius; if d < dmin { dmin = d; mc = kk; } }
					let atom_index = aidx[mc];
					let molecule = atoms[atom_index].molecule;
					let pcen = pijk;
					let outnml = if rp <= 0.0 { point - atoms[atom_index].coor } else { (pcen - point) / rp };
					let other_mol = if molecule == 0 { 1 } else { 0 };
					let mut buried = false;
					for b in atoms.iter() {
						if b.molecule != other_mol { continue; }
						let erl = b.radius + rp;
						let d = pcen.distance_squared(b.coor);
						if d <= erl*erl { buried = true; break; }
					}
					let dot = Dot { coor: point, outnml, area, buried, kind: DotKind::Cavity, atom_index };
					if molecule == 0 { d0.push(dot); } else { d1.push(dot); }
				}
			}
			let n = d0.len() + d1.len();
			if n == 0 { None } else { Some((d0, d1, n)) }
		}).collect();
		for (mut d0, mut d1, n) in results.into_iter() {
			self.run.results.dots.concave += n;
			self.run.dots[0].append(&mut d0);
			self.run.dots[1].append(&mut d1);
		}
		Ok(())
	}
	fn check_probe_collision_idx(&self, point: Vec3, nears: &Vec<usize>, r2: ScValue) -> bool {
		for &np in nears { let p = &self.run.probes[np]; if point.distance_squared(p.point) < r2 { return true; } }
		false
	}

	fn add_dot(&mut self, molecule: usize, kind: DotKind, coor: Vec3, area: ScValue, pcen: Vec3, atom_index: usize) {
		let atom = &self.run.atoms[atom_index];
		let outnml = if self.settings.rp <= 0.0 { coor - atom.coor } else { (pcen - coor) / self.settings.rp };
		let mut buried = false;
		// Robust burial: check against all atoms in the opposite molecule
		let other_mol = if molecule == 0 { 1 } else { 0 };
		for b in self.run.atoms.iter() {
			if b.molecule != other_mol { continue; }
			let erl = b.radius + self.settings.rp;
			let d = pcen.distance_squared(b.coor);
			if d <= erl*erl { buried = true; break; }
		}
		let dot = Dot { coor, outnml, area, buried, kind, atom_index };
		self.run.dots[molecule].push(dot);
	}

	fn distance_point_to_line(&self, cen: Vec3, axis: Vec3, pnt: Vec3) -> ScValue { let vec = pnt - cen; let dt = vec.dot(axis); let mut d2 = vec.magnitude_squared() - dt*dt; if d2 < 0.0 { d2 = 0.0; } d2.sqrt() }

	fn sample_arc(&self, cen: Vec3, rad: ScValue, axis: Vec3, density: ScValue, x: Vec3, v: Vec3, points: &mut Vec<Vec3>) -> Result<ScValue, SurfaceCalculatorError> {
		let y = axis.cross(x);
		let dt1 = v.dot(x);
		let dt2 = v.dot(y);
		let mut angle = dt2.atan2(dt1);
		if angle < 0.0 { angle += 2.0 * PI; }
		self.sample_arc_segment(cen, rad, x, y, angle, density, points)
	}

	fn sample_circle(&self, cen: Vec3, rad: ScValue, axis: Vec3, density: ScValue, points: &mut Vec<Vec3>) -> Result<ScValue, SurfaceCalculatorError> {
		let mut v1 = Vec3::new(axis.y*axis.y + axis.z*axis.z, axis.x*axis.x + axis.z*axis.z, axis.x*axis.x + axis.y*axis.y);
		v1.normalize();
		let dt = v1.dot(axis);
		if dt.abs() > 0.99 { v1 = Vec3::new(1.0, 0.0, 0.0); }
		let mut v2 = axis.cross(v1); v2.normalize();
		let mut x = axis.cross(v2); x.normalize();
		let y = axis.cross(x);
		self.sample_arc_segment(cen, rad, x, y, 2.0*PI, density, points)
	}

	fn sample_arc_segment(&self, cen: Vec3, rad: ScValue, x: Vec3, y: Vec3, angle: ScValue, density: ScValue, points: &mut Vec<Vec3>) -> Result<ScValue, SurfaceCalculatorError> {
		// Match original spacing: delta = 1/(sqrt(density)*rad); sample at midpoints
		if rad <= 0.0 { points.clear(); return Ok(0.0); }
		let delta = 1.0 / (density.sqrt() * rad);
		let mut a = -delta / 2.0;
		points.clear();
		for _ in 0..100000 {
			a += delta;
			if a > angle { break; }
			let c = rad * a.cos();
			let s = rad * a.sin();
			points.push(cen + x*c + y*s);
		}
		if a + delta < angle { return Err(SurfaceCalculatorError::TooManySubdivisions); }
		let ps = if !points.is_empty() { rad * angle / (points.len() as f64) } else { 0.0 };
		Ok(ps)
	}

	pub fn results(&self) -> &Results { &self.run.results }
	pub fn dots(&self, molecule: usize) -> &Vec<Dot> { &self.run.dots[molecule] }

	// Compatibility wrappers (legacy names → new terminology). Safe to remove once callers are updated.
	fn check_point_collision(&self, pcen: Vec3, atoms: &Vec<usize>) -> bool {
		for &idx in atoms.iter().skip(1) {
			let a = &self.run.atoms[idx];
			if pcen.distance(a.coor) <= (a.radius + self.settings.rp) { return true; }
		}
		false
	}
}

// Pure geometry helpers for use in parallel closures (no &self access)
fn geom_sample_arc_segment(cen: Vec3, rad: ScValue, x: Vec3, y: Vec3, angle: ScValue, density: ScValue, points: &mut Vec<Vec3>) -> Result<ScValue, SurfaceCalculatorError> {
	// Match original spacing: delta = 1/(sqrt(density)*rad); sample at midpoints
	if rad <= 0.0 { points.clear(); return Ok(0.0); }
	let delta = 1.0 / (density.sqrt() * rad);
	let mut a = -delta / 2.0;
	points.clear();
	for _ in 0..100000 {
		a += delta;
		if a > angle { break; }
		let c = rad * a.cos();
		let s = rad * a.sin();
		points.push(cen + x*c + y*s);
	}
	if a + delta < angle { return Err(SurfaceCalculatorError::TooManySubdivisions); }
	let ps = if !points.is_empty() { rad * angle / (points.len() as f64) } else { 0.0 };
	Ok(ps)
}

fn geom_sample_arc(cen: Vec3, rad: ScValue, axis: Vec3, density: ScValue, x: Vec3, v: Vec3, points: &mut Vec<Vec3>) -> Result<ScValue, SurfaceCalculatorError> {
	let y = axis.cross(x);
	let dt1 = v.dot(x);
	let dt2 = v.dot(y);
	let mut angle = dt2.atan2(dt1);
	if angle < 0.0 { angle += 2.0 * PI; }
	geom_sample_arc_segment(cen, rad, x, y, angle, density, points)
}

fn geom_sample_circle(cen: Vec3, rad: ScValue, axis: Vec3, density: ScValue, points: &mut Vec<Vec3>) -> Result<ScValue, SurfaceCalculatorError> {
	let mut v1 = Vec3::new(axis.y*axis.y + axis.z*axis.z, axis.x*axis.x + axis.z*axis.z, axis.x*axis.x + axis.y*axis.y);
	v1.normalize();
	let dt = v1.dot(axis);
	if dt.abs() > 0.99 { v1 = Vec3::new(1.0, 0.0, 0.0); }
	let mut v2 = axis.cross(v1); v2.normalize();
	let mut x = axis.cross(v2); x.normalize();
	let y = axis.cross(x);
	geom_sample_arc_segment(cen, rad, x, y, 2.0*PI, density, points)
}
