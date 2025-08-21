use std::ops::{Add, AddAssign, Sub, SubAssign, Mul, Div};

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Vec3 {
	pub x: f64,
	pub y: f64,
	pub z: f64,
}

impl Vec3 {
	pub fn new(x: f64, y: f64, z: f64) -> Self { Self { x, y, z } }
	pub fn zero() -> Self { Self { x: 0.0, y: 0.0, z: 0.0 } }
	pub fn dot(&self, other: Vec3) -> f64 { self.x*other.x + self.y*other.y + self.z*other.z }
	pub fn cross(&self, other: Vec3) -> Vec3 {
		Vec3::new(
			self.y*other.z - self.z*other.y,
			self.z*other.x - self.x*other.z,
			self.x*other.y - self.y*other.x,
		)
	}
	pub fn magnitude_squared(&self) -> f64 { self.dot(*self) }
	pub fn magnitude(&self) -> f64 { self.magnitude_squared().max(0.0).sqrt() }
	pub fn normalize(&mut self) {
		let m = self.magnitude();
		if m > 0.0 { self.x /= m; self.y /= m; self.z /= m; }
	}
	pub fn normalized(&self) -> Vec3 { let mut v = *self; v.normalize(); v }
	pub fn distance_squared(&self, other: Vec3) -> f64 { (*self - other).magnitude_squared() }
	pub fn distance(&self, other: Vec3) -> f64 { self.distance_squared(other).sqrt() }
}

impl Add for Vec3 {
	type Output = Vec3;
	fn add(self, rhs: Vec3) -> Vec3 { Vec3::new(self.x+rhs.x, self.y+rhs.y, self.z+rhs.z) }
}
impl AddAssign for Vec3 {
	fn add_assign(&mut self, rhs: Vec3) { self.x+=rhs.x; self.y+=rhs.y; self.z+=rhs.z; }
}
impl Sub for Vec3 {
	type Output = Vec3;
	fn sub(self, rhs: Vec3) -> Vec3 { Vec3::new(self.x-rhs.x, self.y-rhs.y, self.z-rhs.z) }
}
impl SubAssign for Vec3 {
	fn sub_assign(&mut self, rhs: Vec3) { self.x-=rhs.x; self.y-=rhs.y; self.z-=rhs.z; }
}
impl Mul<f64> for Vec3 {
	type Output = Vec3;
	fn mul(self, rhs: f64) -> Vec3 { Vec3::new(self.x*rhs, self.y*rhs, self.z*rhs) }
}
impl Div<f64> for Vec3 {
	type Output = Vec3;
	fn div(self, rhs: f64) -> Vec3 { Vec3::new(self.x/rhs, self.y/rhs, self.z/rhs) }
}
