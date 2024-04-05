use std::collections::HashMap;
use super::vtk::Mesh;
use super::vtk::Point;
use super::vtk::laplacian_smoothing;
use super::vtk::laplacian_smoothing_with_axis_normalizing;
use super::vtk::laplacian_smoothing_with_center_normalizing;
use super::vtk::laplacian_smoothing_with_cone_normalizing;
use nalgebra as na;
use na::Vector3;

pub struct PocketParameter {
    pub center: Vector3<f32>,
    pub radius: f32,
    pub iteration: i64,
}
pub struct CageParameter {
    pub center: Vector3<f32>,
    pub radius: f32,
    pub axis: Vector3<f32>,
    pub ratio: f32,
    pub theta0: f32,
    pub theta1: f32,
    pub point0: Vector3<f32>,
    pub point1: Vector3<f32>,
    pub iteration: i64,
}
pub struct Cage {
    mesh: Mesh,
    edge_indices: HashMap<String, Vec<i64>>,
    face_indices: HashMap<String, Vec<i64>>,
}
impl Cage {
    pub fn new(mesh: Mesh) -> Self {        
        Cage {
            mesh,
            edge_indices: HashMap::new(),
            face_indices: HashMap::new(),
        }
    }
    pub fn add_edge(&mut self, name: &str, index: Vec<i64>) {
        self.edge_indices.insert(name.to_string(), index);
    }
    pub fn add_face(&mut self, name: &str, index: Vec<i64>) {
        self.face_indices.insert(name.to_string(), index);
    }
    pub fn linspace_curve(&mut self, name: &str, center:Vector3<f32>, radius: f32, theta0: f32, theta1: f32) {

        let curve_index = self.edge_indices.get(name).unwrap();
        let dtheta = (theta1 - theta0) / (curve_index.len() as f32 - 1.0);

        for i in 0..curve_index.len() {
            let theta = theta0 + dtheta * i as f32;

            let x = center[0] + radius * theta.cos();
            let y = center[1] + radius * theta.sin();
            let z = center[2];

            self.mesh.points[curve_index[i] as usize] = Point::new(Vector3::new(x, y, z));
        }
    }
    pub fn linspace_straight(&mut self, name: &str, point0: Vector3<f32>, point1: Vector3<f32>) {
        
        let straight_index = self.edge_indices.get(name).unwrap();

        let dx: Vector3<f32> = (point1 - point0) / (straight_index.len() as f32 - 1.0);

        for i in 0..straight_index.len() {

            let x: Vector3<f32> = point0 + dx * i as f32;

            self.mesh.points[straight_index[i] as usize] = Point::new(x);
        }
    }

    pub fn smooth_curvature(&mut self, name: &str, center:Vector3<f32>, axis: Vector3<f32>, radius: f32, iteration: i64) {
        
        let points = self.mesh.points.clone();
        let curveature_index = self.face_indices.get(name).unwrap();
        let outer_map = self.mesh.outer_map.clone();

        let smoothed_points = laplacian_smoothing_with_axis_normalizing(points, curveature_index.clone(), outer_map, center, axis, radius, iteration);

        for i in curveature_index.clone() {
            self.mesh.points[i as usize] = smoothed_points[i as usize];
        }
    }
    pub fn smooth_sphire(&mut self, name: &str, center:Vector3<f32>, radius: f32, iteration: i64) {
        
        let points = self.mesh.points.clone();
        let sphire_index = self.face_indices.get(name).unwrap();
        let outer_map = self.mesh.outer_map.clone();

        let smoothed_points = laplacian_smoothing_with_center_normalizing(points, sphire_index.clone(), outer_map, center, radius, iteration);

        for i in sphire_index.clone() {
            self.mesh.points[i as usize] = smoothed_points[i as usize];
        }
    }
    pub fn smooth_plane(&mut self, name: &str, iteration: i64) {
        
        let points = self.mesh.points.clone();
        let plane_index = self.face_indices.get(name).unwrap();
        let outer_map = self.mesh.outer_map.clone();

        let smoothed_points = laplacian_smoothing(points, plane_index.clone(), outer_map, iteration);

        for i in plane_index.clone() {
            self.mesh.points[i as usize] = smoothed_points[i as usize];
        }
    }
    pub fn smooth_cone(&mut self, name: &str, center:Vector3<f32>, axis: Vector3<f32>, ratio: f32, iteration: i64) {
        
        let points = self.mesh.points.clone();
        let cone_index = self.face_indices.get(name).unwrap();
        let outer_map = self.mesh.outer_map.clone();

        let smoothed_points = laplacian_smoothing_with_cone_normalizing(points, cone_index.clone(), outer_map, center, axis, ratio, iteration);

        for i in cone_index.clone() {
            self.mesh.points[i as usize] = smoothed_points[i as usize];
        }
    }
    pub fn smooth_inner(&mut self, iteration: i64) {
        
        let points = self.mesh.points.clone();
        let inner_index = self.mesh.inner_index.clone();
        let neighbor_map = self.mesh.neighbor_map.clone();

        let smoothed_points = laplacian_smoothing(points, inner_index.clone(), neighbor_map, iteration);

        for i in inner_index.clone() {
            self.mesh.points[i as usize] = smoothed_points[i as usize];
        }
    }

}

#[cfg(test)]
mod tests {
    // use super::*;
    #[test]
    fn test_lineedge() {
    }
}

