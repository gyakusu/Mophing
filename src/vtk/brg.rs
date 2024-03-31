use std::collections::HashMap;
use super::vtk::laplacian_smoothing_with_cone_normalizing;
use super::vtk::Mesh;
use super::vtk::Point;
use super::vtk::laplacian_smoothing;
use super::vtk::laplacian_smoothing_with_axis_normalizing;
use super::vtk::laplacian_smoothing_with_center_normalizing;

pub struct Cage {
    mesh: Mesh,
    // ellipse_in_center_index: Vec<i64>,
    // ellipse_in_left_index: Vec<i64>,
    // ellipse_in_right_index: Vec<i64>,
    // ellipse_out_center_index: Vec<i64>,
    // ellipse_out_left_index: Vec<i64>,
    // ellipse_out_right_index: Vec<i64>,
    // curve_in_bottom_index: Vec<i64>,
    // curve_in_left_index: Vec<i64>,
    // curve_in_right_index: Vec<i64>,
    // curve_out_bottom_index: Vec<i64>,
    // curve_out_left_index: Vec<i64>,
    // curve_out_right_index: Vec<i64>,
    // curve_bottom_in_index: Vec<i64>,
    // curve_bottom_out_index: Vec<i64>,
    // straight_left_in_index: Vec<i64>,
    // straight_left_out_index: Vec<i64>,
    // straight_right_in_index: Vec<i64>,
    // straight_right_out_index: Vec<i64>,
    // straight_bottom_in_index: Vec<i64>,
    // straight_bottom_out_index: Vec<i64>,
    // straight_indices: HashMap<String, Vec<i64>>,
    // curve_indices: HashMap<String, Vec<i64>>,
    // ellipse_indices: HashMap<String, Vec<i64>>,
    // plane_indices: HashMap<String, Vec<i64>>,
    // curveature_indices: HashMap<String, Vec<i64>>,
    // sphire_indices: HashMap<String, Vec<i64>>,
    // cone_indices: HashMap<String, Vec<i64>>,
    indices: HashMap<String, Vec<i64>>,
}

impl Cage {
    pub fn new(mesh: Mesh) -> Self {        
        Cage {
            mesh: mesh,
            indices: HashMap::new(),
        }
    }
    pub fn add_index(&mut self, name: &str, index: Vec<i64>) {
        self.indices.insert(name.to_string(), index);
    }
    pub fn linspace_curve(&mut self, name: &str, center:[f32; 3], radius: f32, theta0: f32, theta1: f32) {

        let indices = self.indices.get(name).unwrap();
        let dtheta = (theta1 - theta0) / (indices.len() as f32 - 1.0);

        for i in 0..indices.len() {
            let theta = theta0 + dtheta * i as f32;

            let x = center[0] + radius * theta.cos();
            let y = center[1] + radius * theta.sin();
            let z = center[2];

            self.mesh.points[indices[i] as usize] = Point::new([x, y, z]);
        }
    }
    pub fn linspace_straight(&mut self, name: &str, point0: [f32; 3], point1: [f32; 3]) {
        let indices = self.indices.get(name).unwrap();

        let dx: [f32; 3] = (0..3).map(|i| (point1[i] - point0[i]) / (indices.len() as f32 - 1.0)).collect::<Vec<f32>>().try_into().unwrap();

        for i in 0..indices.len() {

            let x: [f32; 3] = (0..3).map(|j| point0[j] + dx[j] * i as f32).collect::<Vec<f32>>().try_into().unwrap(); 

            self.mesh.points[indices[i] as usize] = Point::new(x);
        }
    }
    
    pub fn smooth_curvature(&mut self, name: &str, center:[f32; 3], axis: [f32; 3], radius: f32, iteration: i64) {
        
        let points = self.mesh.points.clone();
        let indices = self.indices.get(name).unwrap();
        let outer_map = self.mesh.outer_map.clone();

        let smoothed_points = laplacian_smoothing_with_axis_normalizing(points, indices.clone(), outer_map, center, axis, radius, iteration);

        for i in indices.clone() {
            self.mesh.points[i as usize] = smoothed_points[i as usize];
        }
    }
    pub fn smooth_sphire(&mut self, name: &str, center:[f32; 3], radius: f32, iteration: i64) {
        
        let points = self.mesh.points.clone();
        let inner_index = self.indices.get(name).unwrap();
        let outer_map = self.mesh.outer_map.clone();

        let smoothed_points = laplacian_smoothing_with_center_normalizing(points, inner_index.clone(), outer_map, center, radius, iteration);

        for i in inner_index.clone() {
            self.mesh.points[i as usize] = smoothed_points[i as usize];
        }
    }
    pub fn smooth_plane(&mut self, name: &str, iteration: i64) {
        
        let points = self.mesh.points.clone();
        let inner_index = self.indices.get(name).unwrap();
        let outer_map = self.mesh.outer_map.clone();

        let smoothed_points = laplacian_smoothing(points, inner_index.clone(), outer_map, iteration);

        for i in inner_index.clone() {
            self.mesh.points[i as usize] = smoothed_points[i as usize];
        }
    }
    pub fn smooth_cone(&mut self, name: &str, center:[f32; 3], axis: [f32; 3], ratio: f32, iteration: i64) {
        
        let points = self.mesh.points.clone();
        let inner_index = self.indices.get(name).unwrap();
        let outer_map = self.mesh.outer_map.clone();

        let smoothed_points = laplacian_smoothing_with_cone_normalizing(points, inner_index.clone(), outer_map, center, axis, ratio, iteration);

        for i in inner_index.clone() {
            self.mesh.points[i as usize] = smoothed_points[i as usize];
        }
    }
    pub fn smooth_inner(&mut self, iteration: i64) {
        
        let points = self.mesh.points.clone();
        let inner_index = self.indices.get("inner").unwrap();
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



// pub trait Edge {
//     fn new(index_size: i64) -> Self;
//     fn as_vec(&self) -> Vec<i64>;
// }
// pub struct LineEdge {
//     index: Vec<i64>,
// }
// pub struct CurveEdge {
//     index: Vec<i64>,
// }
// pub struct PocketEdge {
//     index: Vec<i64>,
// }
// impl LineEdge {
//     pub fn new(index_size: i64) -> Self {
//         LineEdge {
//             index: Vec::with_capacity(index_size),
//         }
//     }
// }
// impl Edge for LineEdge {
//     fn new(index: Vec<i64>) -> Self {
//         LineEdge {
//             index: index,
//         }
//     }
//     fn as_vec(&self) -> Vec<i64> {
//         self.index.clone()
//     }
// }
// impl Edge for CurveEdge {
//     fn new(index: Vec<i64>) -> Self {
//         CurveEdge {
//             index: index,
//         }
//     }
//     fn as_vec(&self) -> Vec<i64> {
//         self.index.clone()
//     }
// }
// impl Edge for PocketEdge {
//     fn new(index: Vec<i64>) -> Self {
//         PocketEdge {
//             index: index,
//         }
//     }
//     fn as_vec(&self) -> Vec<i64> {
//         self.index.clone()
//     }
// }
// pub struct CurveFace {
// }

// pub struct SphireFace {
// }

// pub trait Face {
//     fn new() -> Self;
// }
    // pub fn add_curve(&mut self, name: &str, indices: Vec<i64>) {
    //     self.indices.insert(name.to_string(), indices);
    // }
    // pub fn add_straight(&mut self, name: &str, indices: Vec<i64>) {
    //     self.indices.insert(name.to_string(), indices);
    // }
    // pub fn add_plane(&mut self, name: &str, indices: Vec<i64>) {
    //     self.plane_indices.insert(name.to_string(), indices);
    // }
    // pub fn add_curveature(&mut self, name: &str, indices: Vec<i64>) {
    //     self.curveature_indices.insert(name.to_string(), indices);
    // }
    // pub fn add_sphire(&mut self, name: &str, indices: Vec<i64>) {
    //     self.sphire_indices.insert(name.to_string(), indices);
    // }
    // pub fn add_cone(&mut self, name: &str, indices: Vec<i64>) {
    //     self.cone_indices.insert(name.to_string(), indices);
    // }
