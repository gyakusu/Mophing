use std::f32::consts::PI;
use std::vec;
use std::collections::HashMap;
use std::collections::HashSet;
use pyo3::prelude::*;

use nalgebra as na;
use na::Vector3;

use super::point::Point;
use super::mesh::Mesh;
use super::io::read_vtk;
use super::io::read_edge_from_xml;
use super::io::read_index_from_xml;
use super::io::copy_vtk_and_replace_point;

use super::mesh::laplacian_smoothing_on_plane;
use super::mesh::laplacian_smoothing_with_axis_normalizing;
use super::mesh::laplacian_smoothing_with_center_normalizing;
use super::mesh::laplacian_smoothing_with_cotangent_and_center_normalizing;
use super::mesh::laplacian_smoothing_wrapper;
// use super::mesh::laplacian_smoothing_with_cone_normalizing;

fn linspace_arc(i: usize, center: Vector3<f32>, radius: f32, theta: f32, dtheta: f32) -> Point {
    let t: f32 = theta + dtheta * (i as f32);
    Point::new(
        center[0] + radius * t.sin(),
        center[1] + radius * t.cos(),
        center[2]
    )
}
fn linspace_line(i: usize, x: Vector3<f32>, dx: Vector3<f32>) -> Point {
        Point::from_vec(x + dx * (i as f32))
}
// 円柱周りの（近似）楕円軌道を生成する関数．
fn linspace_ellipse(i: usize, x: f32, z: f32, r: f32, dx: f32, dz: f32, theta: f32, dtheta: f32) -> Point {
    let t = theta + dtheta * (i as f32);

    let x0 = x + t.sin() * dx;
    let y0 = (r * r - x0 * x0).sqrt();
    let z0 = z + t.cos() * dz;

    Point::new(x0, y0, z0)
}
fn row_of_cosine(a: f32,b: f32,c: f32) -> f32 {
    (b*b + c*c - a*a) / (2.0 * b * c)
}
fn calculate_x_diameter(r0: f32, y: f32, r1: f32) -> (f32, f32) {
    let cos_r1 = row_of_cosine(r1, r0, y);
    let sin_r1 = (1.0 - cos_r1 * cos_r1).sqrt();
    let dx = r0 * sin_r1;
    let dy = r0 * cos_r1;
    (dx, dy)
}
fn calculate_z_diameter(r0: f32, y: f32, r1: f32) -> f32 {
    let dy = y - r0;
    let dz = (r1 * r1 - dy * dy).sqrt();
    dz
}
fn cutted_sphire_radius(r: f32, h: f32) -> f32 {
    (r*r - h*h).sqrt()
}
// メッシュ構造体パラメタ
#[derive(Clone, Debug)]
pub struct PocketParameter {
    pub x: Vector3<f32>,
    pub r: f32,
}
#[derive(Clone, Debug)]
pub struct NeckParameter {
    pub x: Vector3<f32>,
    pub r: f32,
    pub h: f32,
    pub dh: f32,      // 保持器最高点からポケット面までの距離
    pub h_ratio: f32, // 面取までの距離／dhの比率
    pub r_ratio: f32, // 面取半径／dhの比率
}
#[pyclass]
#[derive(Clone, Debug)]
pub struct CageParameter {
    pub axis: Vector3<f32>,
    pub theta0: f32,
    pub theta1: f32,
    pub r0: f32,
    pub r1: f32,
    pub h0: f32,
    pub h1: f32,
    pub bevel: f32,

    pub pocket: PocketParameter,
    pub neck: NeckParameter,
    // pub face_list: Vec<String>,
}
impl CageParameter {
    pub fn sample() -> Self {
        CageParameter {
            axis: Vector3::new(0.0, 0.0, 1.0),
            theta0: -PI / 6.0,
            theta1:  PI / 6.0,
            r0: 2.345e-3,
            r1: 2.850e-3,
            h0: 0.93e-3,
            h1: 2.10e-3,
            bevel: 0.10e-3,
            
            pocket: PocketParameter {
                x: Vector3::new(0.0, 2.65e-3, 2.00e-3),
                r: 0.825e-3,
            },
            neck: NeckParameter {
                x: Vector3::new(0.0, 2.65e-3, 2.00e-3 - 0.30e-3),
                r: 1.20e-3,
                h: 2.45e-3,
                dh: 0.1519e-3,
                h_ratio: 0.355,
                r_ratio: 0.55,
            },
        }
    }
}
#[pymethods]
impl CageParameter {
    #[new]
    pub fn new(r0: f32, r1: f32, h0: f32, h1: f32, bevel: f32, pr: f32, nz: f32, nr: f32, nh: f32, ndh: f32) -> Self {
        let pocket = PocketParameter {
            x: Vector3::new(0.0, 2.65e-3, 2.00e-3),
            r: pr,
        };
        let neck = NeckParameter {
            x: Vector3::new(0.0, 2.65e-3, nz),
            r: nr,
            h: nh,
            dh: ndh,
            h_ratio: 0.355,
            r_ratio: 0.55,
        };
        let cage_parameter = CageParameter {
            axis: Vector3::new(0.0, 0.0, 1.0),
            theta0: -PI / 6.0,
            theta1:  PI / 6.0,
            r0,
            r1,
            h0,
            h1,
            bevel,
            pocket,
            neck,
        };
        cage_parameter
    }
}
#[derive(Clone, Debug)]
pub struct BallParameter {
    pub r: f32,
    pub x: Vector3<f32>,
}
impl BallParameter {
    pub fn sample() -> Self {
        BallParameter {
            r: 0.76e-3,
            // r: 0.778e-3,
            x: Vector3::new(0.0, 2.645e-3, 2.0e-3),
            // map: HashMap::new(),
        }
    }
}
#[derive(Clone, Debug)]
#[pyclass]
pub struct Brg {
    mesh: Mesh,
    cage: CageParameter,
    ball: BallParameter,

    edges: HashMap<String, Vec<usize>>,
    faces: HashMap<String, HashSet<usize>>,
    periodics: (HashSet<usize>, HashSet<usize>),
    
}
#[pymethods]
impl Brg {
    #[new]
    pub fn from_file(vtk_path: &str, index_path: &str, cage: &CageParameter) -> Self {        
        let mesh = read_vtk(vtk_path);
        let edges = read_edge_from_xml(index_path, "edge").unwrap();
        let faces = read_index_from_xml(index_path, "face").unwrap();
        let mut brg = Brg::new(&mesh, cage);
        brg.set_edge_and_face(&edges, &faces);
        brg
    }
    pub fn reload_cage(&mut self, cage: &CageParameter, i: usize, j: usize) {
        self.cage = cage.clone();
        self.linspace_all();
        for _ in 0..i {
            self.smooth_face();
        }
        for _ in 0..j {
            self.smooth_ball();
            self.smooth_inner();
        }
    }
    pub fn get_points_as_list(&self, py: Python) -> PyResult<PyObject> {
        let points = self.get_points(); // Assuming this returns Vec<[f64; 3]>
        let py_points: Vec<Vec<f64>> = points.iter().map(|p| p.as_vec().iter().map(|&x| x as f64).collect()).collect();
        Ok(py_points.into_py(py))
    }
    pub fn extract_surface_as_list(&self, py: Python) -> PyResult<(PyObject, PyObject)> {
        let (reindex_map, surface_faces) = self.mesh.extract_surface();

        let mut py_reindex_map: Vec<usize> = vec![0; reindex_map.len()];
        for (old, new) in reindex_map.iter() {
            py_reindex_map[*new] = *old;
        }
        let mut py_faces: Vec<Vec<usize>>  = vec![vec![0; 3]; surface_faces.len()];
        for (i, (face, is_front)) in surface_faces.iter().enumerate() {
            let face_index = face.vec_with_order(*is_front);
            for j in 0..3 {
                py_faces[i][j] = *reindex_map.get(&face_index[j]).unwrap();
            }
        }
        Ok((py_reindex_map.into_py(py), py_faces.into_py(py)))
    }
    pub fn smooth_face(&mut self) {

        let points = self.get_points();
        let bottom: Vector3<f32>    = Vector3::new(0.0, 0.0, self.cage.h0);
        let shoulder: Vector3<f32>  = Vector3::new(0.0, 0.0, self.cage.h1);
        let top: Vector3<f32>       = Vector3::new(0.0, 0.0, self.cage.neck.h);
        let axis: Vector3<f32>      = Vector3::new(0.0, 0.0, 1.0);
        let axis_left: Vector3<f32> = Vector3::new(self.cage.theta0.cos(), self.cage.theta0.sin(), 0.0);
        let axis_right: Vector3<f32> = Vector3::new(self.cage.theta1.cos(), self.cage.theta1.sin(), 0.0);

        let get_inner_index = |name: &str| self.faces.get(name).unwrap_or(&HashSet::new()).clone();
        let surface_map = self.mesh.surface_map.clone();

        let operations: [(&str, fn(&Vec<Point>, &HashSet<usize>, &HashMap<usize, HashSet<usize>>, Vector3<f32>, Vector3<f32>, f32) -> Vec<Point>, Vector3<f32>, Vector3<f32>, f32); 18] =
        [
            ("curvature_in", laplacian_smoothing_with_axis_normalizing 
            as fn(&_, &_, &_, _, _, _) -> _, bottom, axis, self.cage.r0),

            ("curvature_out", laplacian_smoothing_with_axis_normalizing 
            as fn(&_, &_, &_, _, _, _) -> _, bottom, axis, self.cage.r1),

            ("bottom", laplacian_smoothing_on_plane 
            as fn(&_, &_, &_, _, _, _) -> _, bottom, axis, 0.),

            ("top_left", laplacian_smoothing_on_plane 
            as fn(&_, &_, &_, _, _, _) -> _, top, axis, 0.),

            ("top_right", laplacian_smoothing_on_plane 
            as fn(&_, &_, &_, _, _, _) -> _, top, axis, 0.),

            ("shoulder_left", laplacian_smoothing_on_plane 
            as fn(&_, &_, &_, _, _, _) -> _, shoulder, axis, 0.),

            ("shoulder_right", laplacian_smoothing_on_plane 
            as fn(&_, &_, &_, _, _, _) -> _, shoulder, axis, 0.),

            ("sphire", laplacian_smoothing_with_center_normalizing 
            as fn(&_, &_, &_, _, _, _) -> _, self.cage.pocket.x, axis, self.cage.pocket.r),

            ("sphire_left", laplacian_smoothing_with_center_normalizing 
            as fn(&_, &_, &_, _, _, _) -> _, self.cage.neck.x, axis, self.cage.neck.r),

            ("sphire_right", laplacian_smoothing_with_center_normalizing 
            as fn(&_, &_, &_, _, _, _) -> _, self.cage.neck.x, axis, self.cage.neck.r),

            ("periodic_left", laplacian_smoothing_on_plane 
            as fn(&_, &_, &_, _, _, _) -> _, bottom, axis_left, 0.),

            ("periodic_right", laplacian_smoothing_on_plane 
            as fn(&_, &_, &_, _, _, _) -> _, bottom, axis_right, 0.),

            ("apature_left", laplacian_smoothing_wrapper
            as fn(&_, &_, &_, _, _, _) -> _, axis, axis, 0.),

            ("apature_right", laplacian_smoothing_wrapper 
            as fn(&_, &_, &_, _, _, _) -> _, axis, axis, 0.),

            ("cone_in", laplacian_smoothing_wrapper 
            as fn(&_, &_, &_, _, _, _) -> _, axis, axis, 0.),

            ("cone_out", laplacian_smoothing_wrapper 
            as fn(&_, &_, &_, _, _, _) -> _, axis, axis, 0.),

            ("straw_left", laplacian_smoothing_wrapper 
            as fn(&_, &_, &_, _, _, _) -> _, axis, axis, 0.),

            ("straw_right", laplacian_smoothing_wrapper 
            as fn(&_, &_, &_, _, _, _) -> _, axis, axis, 0.),            
        ];
        for (name, function, arg1, arg2, arg3) in operations.iter() {
            let inner_index = get_inner_index(name);
            let smoothed_point = function(&points, &inner_index, &surface_map, *arg1, *arg2, *arg3);
            inner_index.iter().for_each(|&i| self.mesh.points[i] = smoothed_point[i]);
        }
    }
    pub fn smooth_ball(&mut self) {
        let inner_index = self.faces.get("on_ball").unwrap_or(&HashSet::new()).clone();
        let neighbor_map = self.mesh.neighbor_map.clone();
        let new_points = laplacian_smoothing_with_center_normalizing(&self.get_points(), &inner_index, &neighbor_map, self.ball.x, Vector3::zeros(), self.ball.r);
        for i in inner_index.clone() {
            self.mesh.points[i] = new_points[i];
        }
    }
    pub fn smooth_ball_with_cotangent(&mut self) {
        let inner_index = self.faces.get("on_ball").unwrap_or(&HashSet::new()).clone();
        let neighbor_map = self.mesh.neighbor_map.clone();
        let new_points = laplacian_smoothing_with_cotangent_and_center_normalizing(&self.get_points(), &inner_index, &neighbor_map, self.ball.x, Vector3::zeros(), self.ball.r);
        for i in inner_index.clone() {
            self.mesh.points[i] = new_points[i];
        }
    }
    pub fn smooth_inner(&mut self) {
        self.mesh.smooth_inner();
    }
    pub fn smooth_inner_with_cotangent(&mut self) {
        self.mesh.smooth_inner_with_cotangent();
    }
    pub fn smooth_inner_with_flower(&mut self) {
        self.mesh.smooth_inner_with_flower();
    }
    pub fn linspace_all(&mut self) {
        for name in self.edges.keys() {
            let index = self.edges.get(name).unwrap();
            let points: Vec<Point> = self.linspace(name);

            #[cfg(debug_assertions)] {
                let old_points_all = self.get_points();
                let old_points: Vec<Point> = index.iter().map(|i| old_points_all[*i as usize]).collect();
                let distance_small = (points.first().unwrap().as_vec() - old_points.first().unwrap().as_vec()).norm();
                let distance_large = (points.first().unwrap().as_vec() - old_points.last().unwrap().as_vec()).norm();
                if distance_small > distance_large {
                    // panic!("must be flipped at {}", name)
                }
                // println!("{}", distance_small); // 目視で確認．点の移動が大きい場合はflip_xを使う．
            }
            for i in 0..index.len() {
                self.mesh.points[index[i] as usize] = points[i];
            }
        }
    }
    pub fn scale(&mut self, s: f32) {
        self.mesh.points.iter_mut().for_each(|p| *p = p.mul(s));
    }
    pub fn std(&self) -> f32 {
        let points_vec = self.get_points();
        let mean_of_points = points_vec.iter().fold(Vector3::zeros(), |sum, p| sum + p.as_vec()) / (points_vec.len() as f32);
        let variance_of_points = points_vec.iter().fold(0.0, |sum, p| sum + (p.as_vec() - mean_of_points).norm());
        let s = (variance_of_points / (points_vec.len() as f32)).sqrt();
        s
    }
    pub fn normalize_center(&mut self) {
        let ball_index = self.faces.get("on_ball").unwrap_or(&HashSet::new()).clone();
        self.mesh.normalize_center(self.ball.x, self.ball.r, &ball_index);

        let sphire_index = self.faces.get("sphire").unwrap_or(&HashSet::new()).clone();
        self.mesh.normalize_center(self.cage.pocket.x, self.cage.pocket.r, &sphire_index);
    }
    pub fn flip_negative_volume(&mut self, ratio: f32) -> bool {
        let ball_index: HashSet<usize> = self.faces.get("on_ball").unwrap_or(&HashSet::new()).clone();
        let inner_index: HashSet<usize> = self.mesh.inner_index.union(&ball_index).cloned().collect();

        let sphire_index: HashSet<usize> = self.faces.get("sphire").unwrap_or(&HashSet::new()).clone();
        let inner_index: HashSet<usize> = inner_index.union(&sphire_index).cloned().collect();

        self.mesh.flip_negative_volume(&inner_index, ratio, true)
    }
    pub fn flip_negative_volume_only_sphire(&mut self, ratio: f32) -> bool {
        let ball_index: HashSet<usize> = self.faces.get("on_ball").unwrap_or(&HashSet::new()).clone();
        let sphire_index: HashSet<usize> = self.faces.get("sphire").unwrap_or(&HashSet::new()).clone();
        let inner_index: HashSet<usize> = ball_index.union(&sphire_index).cloned().collect();

        self.mesh.flip_negative_volume(&inner_index, ratio, false)
    }
    pub fn write_vtk_from_base(&self, origin_path: &str, write_path: &str) {
        copy_vtk_and_replace_point(origin_path, write_path, &self.get_points());
    }
    pub fn flip_with_scale(&mut self, scale: f32, ratio: f32) -> bool {
        self.scale(1.0 / scale);
        let result = self.flip_negative_volume(ratio);
        self.scale(scale);
        self.normalize_center();
        result
    }
    pub fn flip_sphire_with_scale(&mut self, scale: f32, ratio: f32) -> bool {
        self.scale(1.0 / scale);
        let result = self.flip_negative_volume_only_sphire(ratio);
        self.scale(scale);
        self.normalize_center();
        result
    }

}
impl Brg {
    pub fn new(mesh: &Mesh, cage: &CageParameter) -> Self {        
        Brg {
            mesh: mesh.clone(),
            cage: cage.clone(),
            ball: BallParameter::sample(),
            edges: HashMap::new(),   // エッジのインデックス
            faces: HashMap::new(),   // エッジを除く面のインデックス
            periodics: (HashSet::new(), HashSet::new()), // 周期境界のインデックス
            
        }
    }    pub fn sample(mesh: &Mesh) -> Self {
        let cage = CageParameter::sample();
        Brg::new(mesh, &cage)
    }
    pub fn set_edge_and_face(&mut self, edges: &HashMap<String, Vec<usize>>, faces: &HashMap<String, HashSet<usize>>) {
        self.edges = edges.clone();
        self.faces = faces.clone();
    }
    pub fn set_periodic(&mut self, left: &HashSet<usize>, right: &HashSet<usize>) {
        self.periodics = (left.clone(), right.clone());
    }
    pub fn get_points(&self) -> Vec<Point> {
        self.mesh.points.clone()
    }
    pub fn set_points(&mut self, points: Vec<Point>) {
        self.mesh.points = points;
    }
    fn get_curve_params(&self, cage: &CageParameter, h: f32, n: usize) -> (Vector3<f32>, f32, f32) {
        let center: Vector3<f32> = Vector3::new(0.0, 0.0, h);
        let theta = cage.theta0;
        let dtheta = (cage.theta1 - cage.theta0) / (n as f32 - 1.0);
        (center, theta, dtheta)
    }
    fn get_asymmetry_curve_params(&self, cage: &CageParameter, h: f32, n: usize, r0: f32, r1: f32) -> (Vector3<f32>, f32, f32, f32) {
        let center: Vector3<f32> = Vector3::new(0.0, 0.0, h);
        let y = cage.neck.x[1];
        let theta = cage.theta0;
        let theta1 = -(calculate_x_diameter(r0, y, r1).0 / r0).asin();
        let dtheta = (theta1 - theta) / (n as f32 - 1.0);
        let radius = r0;
        (center, theta, dtheta, radius)
    }
    fn get_asymmetry_curve_top_params(&self, cage: &CageParameter, h: f32, n: usize, r0: f32, r1: f32, r2: f32) -> (Vector3<f32>, f32, f32, f32) {
        let center: Vector3<f32> = Vector3::new(0.0, 0.0, h);
        let y = cage.pocket.x[1];
        let theta = (calculate_x_diameter(r0, y, r1).0 / r0).asin();
        let theta1 = (calculate_x_diameter(r0, y, r2).0 / r0).asin();
        let dtheta = (theta1 - theta) / (n as f32 - 1.0);
        let radius = r0;
        (center, theta, dtheta, radius)
    }
    fn get_collar_params(&self, cage: &CageParameter, h: f32, n: usize, r: f32, x: f32, y: f32) -> (Vector3<f32>, f32, f32, f32) {
        let center: Vector3<f32> = Vector3::new(x, y, h);
        let y0 = calculate_x_diameter(cage.r0, y, r).1;
        let y1 = calculate_x_diameter(cage.r1, y, r).1;
        let theta = ((y0 - y) / r).acos();
        let theta1 = ((y1 - y) / r).acos();
        let dtheta = (theta1 - theta) / (n as f32 - 1.0);
        let radius = r;
        (center, -theta, dtheta, radius)
    }
    fn get_ellipse_params(&self, cage: &CageParameter, r: f32, n: usize) -> (f32, f32, f32, f32, f32, f32, f32) {
        let x = cage.pocket.x[0];
        let y = cage.pocket.x[1];
        let z = cage.pocket.x[2];
        let dx = calculate_x_diameter(r, y, cage.pocket.r).0;
        let dz = calculate_z_diameter(r, y, cage.pocket.r);
        let cos_theta = (cage.neck.h - cage.neck.dh - z) / dz;
        let theta = cos_theta.acos();
        let dtheta = 2.0 * (PI - theta) / (n as f32 - 1.0);        
        (x, z, r, dx, dz, theta, dtheta)
    }
    pub fn get_settings(&self) -> (HashSet<usize>, HashMap<usize, HashSet<usize>>, HashMap<usize, HashSet<usize>>) {

        let inner_index = self.mesh.inner_index.clone();
        let neighbor_map = self.mesh.neighbor_map.clone();
        let surface_map = self.mesh.surface_map.clone();

        (inner_index, neighbor_map, surface_map)
    }
    pub fn linspace(&self, name: &str) -> Vec<Point> {
    
        let cage = &self.cage;
        let n = self.edges.get(name).unwrap().len();

        let r_middle = cutted_sphire_radius(cage.pocket.r, cage.neck.h - cage.pocket.x[2] - cage.neck.dh); // 円筒部における半径
        let r_neck   = cutted_sphire_radius(cage.neck.r, cage.h1 - cage.neck.x[2]); // 肩部における爪部の半径
        let r_top    = cutted_sphire_radius(cage.neck.r, cage.neck.h - cage.neck.x[2]); // 頂部における爪外径の半径
        let r_apt    = r_middle + cage.neck.dh * cage.neck.r_ratio; // 頂部における爪内径の半径

        // let r_cone = 1.0; // 保持器底面の面取りの傾斜．内径・外径共通にしてる．
        // let x_apt: Vector3<f32> = Vector3::new(cage.pocket.x[0], cage.pocket.x[1], cage.neck.h - cage.neck.dh - cage.pocket.r / cage.neck.h_ratio);
        // let x_cone_in: Vector3<f32>  = Vector3::new(0., 0., cage.h0 + (cage.r0 + cage.bevel) / r_cone);
        // let x_cone_out: Vector3<f32> = Vector3::new(0., 0., cage.h0 - (cage.r1 - cage.bevel) / r_cone);

        let mut points: Vec<Point> = Vec::with_capacity(n);

        match name {
            "curve_bottom_in" => {
                let (center, theta, dtheta) = self.get_curve_params(cage, cage.h0, n);
                (0..n).for_each(|i| points.push(linspace_arc(i, center, cage.r0 + cage.bevel, theta, dtheta)));
            },
            "curve_bottom_out" => {
                let (center, theta, dtheta) = self.get_curve_params(cage, cage.h0, n);
                (0..n).for_each(|i| points.push(linspace_arc(i, center, cage.r1 - cage.bevel, theta, dtheta)));
            },
            "curve_in_bottom" => {
                let (center, theta, dtheta) = self.get_curve_params(cage, cage.h0 + cage.bevel, n);
                (0..n).for_each(|i| points.push(linspace_arc(i, center, cage.r0, theta, dtheta)));
            },
            "curve_out_bottom" => {
                let (center, theta, dtheta) = self.get_curve_params(cage, cage.h0 + cage.bevel, n);
                (0..n).for_each(|i| points.push(linspace_arc(i, center, cage.r1, theta, dtheta)));
            },
            "curve_in_top_left" | "curve_in_top_right" => {
                let (center, theta, dtheta, radius) = self.get_asymmetry_curve_params(cage, cage.h1, n, cage.r0, r_neck);
                (0..n).for_each(|i| points.push(linspace_arc(i, center, radius, theta, dtheta)));
            },
            "curve_out_top_left" | "curve_out_top_right" => {
                let (center, theta, dtheta, radius) = self.get_asymmetry_curve_params(cage, cage.h1, n, cage.r1, r_neck);
                (0..n).for_each(|i| points.push(linspace_arc(i, center, radius, theta, dtheta)));
            },
            "curve_top_in_left" | "curve_top_in_right" => {
                let (center, theta, dtheta, radius) = self.get_asymmetry_curve_top_params(cage, cage.neck.h, n, cage.r0, r_top, r_apt);
                (0..n).for_each(|i| points.push(linspace_arc(i, center, radius, -theta, -dtheta)));
            }
            "curve_top_out_left" | "curve_top_out_right" => {
                let (center, theta, dtheta, radius) = self.get_asymmetry_curve_top_params(cage, cage.neck.h, n, cage.r1, r_top, r_apt);
                (0..n).for_each(|i| points.push(linspace_arc(i, center, radius, -theta, -dtheta)));
            },
            "straight_left_bottom" | "straight_right_bottom" => {
                let r0 = cage.r0 + cage.bevel;
                let dr = (cage.r1 - cage.bevel - r0) / (n as f32 - 1.0);
                let sin = cage.theta0.sin();
                let cos = cage.theta0.cos();
                let x: Vector3<f32> = Vector3::new(r0 * sin, r0 * cos, cage.h0);
                let dx: Vector3<f32> = dr * Vector3::new(sin, cos, 0.0);
                (0..n).for_each(|i| points.push(linspace_line(i, x, dx)));
            },
            "straight_in_left" | "straight_in_right" => {
                let x = cage.r0 * cage.theta0.sin();
                let y = cage.r0 * cage.theta0.cos();
                let x0: Vector3<f32> = Vector3::new(x, y, cage.h0 + cage.bevel);
                let x1: Vector3<f32> = Vector3::new(x, y, cage.h1);
                let dx: Vector3<f32> = (x1 - x0) / (n as f32 - 1.0);
                (0..n).for_each(|i| points.push(linspace_line(i, x0, dx)));
            },
            "straight_out_left" | "straight_out_right" => {
                let x = cage.r1 * cage.theta0.sin();
                let y = cage.r1 * cage.theta0.cos();
                let x0: Vector3<f32> = Vector3::new(x, y, cage.h0 + cage.bevel);
                let x1: Vector3<f32> = Vector3::new(x, y, cage.h1);
                let dx: Vector3<f32> = (x1 - x0) / (n as f32 - 1.0);
                (0..n).for_each(|i| points.push(linspace_line(i, x0, dx)));
            },
            "straight_top_left" | "straight_top_right" => {
                let sin = cage.theta0.sin();
                let cos = cage.theta0.cos();
                let (r0, r1) = (cage.r0, cage.r1);
                let dr = (r1 - r0) / (n as f32 - 1.0);
                let x: Vector3<f32> = Vector3::new(r0 * sin, r0 * cos, cage.h1);
                let dx: Vector3<f32> = dr * Vector3::new(sin, cos, 0.0);
                (0..n).for_each(|i| points.push(linspace_line(i, x, dx)));
            },
            "slope_bottom_left_in" | "slope_bottom_right_in" => {
                let sin = cage.theta0.sin();
                let cos = cage.theta0.cos();
                let (r0, r1) = (cage.r0, cage.r0 + cage.bevel);
                let x0: Vector3<f32> = Vector3::new(r1 * sin, r1 * cos, cage.h0);
                let x1: Vector3<f32> = Vector3::new(r0 * sin, r0 * cos, cage.h0 + cage.bevel);
                let dx: Vector3<f32> = (x1 - x0) / (n as f32 - 1.0);
                (0..n).for_each(|i| points.push(linspace_line(i, x0, dx)));
            }
            "slope_bottom_left_out" | "slope_bottom_right_out" => {
                let sin = cage.theta0.sin();
                let cos = cage.theta0.cos();
                let(r0, r1) = (cage.r1 - cage.bevel, cage.r1);
                let x0: Vector3<f32> = Vector3::new(r0 * sin, r0 * cos, cage.h0);
                let x1: Vector3<f32> = Vector3::new(r1 * sin, r1 * cos, cage.h0 + cage.bevel);
                let dx: Vector3<f32> = (x1 - x0) / (n as f32 - 1.0);
                (0..n).for_each(|i| points.push(linspace_line(i, x0, dx)));
            }
            "slope_in_top_left" | "slope_in_top_right" => {
                let y = cage.neck.x[1];
                let (x0, y0) = calculate_x_diameter(cage.r0, y, r_apt);
                let (x1, y1) = calculate_x_diameter(cage.r0, y, r_middle);
                let xyz1: Vector3<f32> = Vector3::new(-x0, y0, cage.neck.h);
                let xyz0: Vector3<f32> = Vector3::new(-x1, y1, cage.neck.h - cage.neck.dh * cage.neck.h_ratio);
                let dx: Vector3<f32> = (xyz1 - xyz0) / (n as f32 - 1.0);
                (0..n).for_each(|i| points.push(linspace_line(i, xyz0, dx)));
            }
            "slope_out_top_left" | "slope_out_top_right" => {
                let y = cage.neck.x[1];
                let (x0, y0) = calculate_x_diameter(cage.r1, y, r_apt);
                let (x1, y1) = calculate_x_diameter(cage.r1, y, r_middle);
                let xyz1: Vector3<f32> = Vector3::new(-x0, y0, cage.neck.h);
                let xyz0: Vector3<f32> = Vector3::new(-x1, y1, cage.neck.h - cage.neck.dh * cage.neck.h_ratio);
                let dx: Vector3<f32> = (xyz1 - xyz0) / (n as f32 - 1.0);
                (0..n).for_each(|i| points.push(linspace_line(i, xyz0, dx)));
            }
            "collar_left_out" | "collar_right_out" => {
                let (center, theta, dtheta, radius) = self.get_collar_params(cage, cage.h1, n, r_neck, cage.neck.x[0], cage.neck.x[1]);
                (0..n).for_each(|i| points.push(linspace_arc(i, center, radius, theta, -dtheta)));
            }
            "collar_left_middle" | "collar_right_middle" => {
                let (center, theta, dtheta, radius) = self.get_collar_params(cage, cage.neck.h, n, r_top, cage.neck.x[0], cage.neck.x[1]);
                (0..n).for_each(|i| points.push(linspace_arc(i, center, radius, theta, -dtheta)));
            }
            "collar_left_in" | "collar_right_in" => {
                let (center, theta, dtheta, radius) = self.get_collar_params(cage, cage.neck.h, n, r_apt, cage.neck.x[0], cage.neck.x[1]);
                (0..n).for_each(|i| points.push(linspace_arc(i, center, radius, theta, -dtheta)));
            },
            "collar_bottom_left" | "collar_bottom_right" => {
                let (center, theta, dtheta, radius) = self.get_collar_params(cage, cage.neck.h - cage.neck.dh, n, r_middle, cage.pocket.x[0], cage.pocket.x[1]);
                (0..n).for_each(|i| points.push(linspace_arc(i, center, radius, theta, -dtheta)));
            },
            "collar_middle_left" | "collar_middle_right" => {
                let (center, theta, dtheta, radius) = self.get_collar_params(cage, cage.neck.h - cage.neck.dh * cage.neck.h_ratio, n, r_middle, cage.pocket.x[0], cage.pocket.x[1]);
                (0..n).for_each(|i| points.push(linspace_arc(i, center, radius, theta, -dtheta)));
            },
            "ellipse_in_center" => {
                let (x, z, r, dx, dz, theta, dtheta) = self.get_ellipse_params(cage, cage.r0, n);
                (0..n).for_each(|i| points.push(linspace_ellipse(i, x, z, r, dx, dz, -theta, -dtheta)));
            },
            "ellipse_out_center" => {
                let (x, z, r, dx, dz, theta, dtheta) = self.get_ellipse_params(cage, cage.r1, n);
                (0..n).for_each(|i| points.push(linspace_ellipse(i, x, z, r, dx, dz, -theta, -dtheta)));
            },
            "ellipse_in_left" | "ellipse_in_right" => {
                let z0 = cage.h1;
                let z1 = cage.neck.h;
                let dz = (z1 - z0) / (n as f32 - 1.0);
                for i in 0..n {
                    let h = z0 + dz * (i as f32);
                    let r = cutted_sphire_radius(self.cage.neck.r, h - cage.neck.x[2]);
                    let (x, y) = calculate_x_diameter(self.cage.r0, self.cage.neck.x[1], r);
                    points.push(Point::new(-x, y, h));
                }
            },
            "ellipse_out_left" | "ellipse_out_right" => {
                let z0 = cage.h1;
                let z1 = cage.neck.h;
                let dz = (z1 - z0) / (n as f32 - 1.0);
                for i in 0..n {
                    let h = z0 + dz * (i as f32);
                    let r = cutted_sphire_radius(self.cage.neck.r, h - cage.neck.x[2]);
                    let (x, y) = calculate_x_diameter(self.cage.r1, self.cage.neck.x[1], r);
                    points.push(Point::new(-x, y, h));
                }
            },
            "stripe_in_left" | "stripe_in_right" => {
                let (x, y) = calculate_x_diameter(cage.r0, cage.pocket.x[1], r_middle);
                let z0 = cage.neck.h - cage.neck.dh;
                let z1 = cage.neck.h - cage.neck.dh * cage.neck.h_ratio;
                let dz = (z1 - z0) / (n as f32 - 1.0);
                (0..n).for_each(|i| points.push(Point::new(-x, y, z0 + dz * (i as f32))));
            },
            "stripe_out_left" | "stripe_out_right" => {
                let (x, y) = calculate_x_diameter(cage.r1, cage.pocket.x[1], r_middle);
                let z0 = cage.neck.h - cage.neck.dh;
                let z1 = cage.neck.h - cage.neck.dh * cage.neck.h_ratio;
                let dz = (z1 - z0) / (n as f32 - 1.0);
                (0..n).for_each(|i| points.push(Point::new(-x, y, z0 + dz * (i as f32))));
            },
            _ => {
            }
        }
        let flip_list = vec![
            "curve_in_top_right", 
            "curve_out_top_right", 
            "slope_in_top_left",
            "slope_out_top_left",
            "curve_top_in_right",
            "curve_top_out_right",
            ];
        if flip_list.contains(&name) {
            points.reverse();
        }
        if name.contains("right") {
            return points.iter().map(|p| p.flip_x()).collect();
        }
        points
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::io;

    #[test]
    fn test_calculate_minor_diameter() {
        let r0 = 4.0;
        let y = 5.0;
        let r1 = 3.0;
        let (dx, _) = calculate_x_diameter(r0, y, r1);
        assert_ne!(dx, r0 * 0.6);
    }
    #[test]
    fn test_calculate_major_diameter() {
        let r0 = 7.0;
        let y = 4.0;
        let r1 = 5.0;

        let dz = calculate_z_diameter(r0, y, r1);
        assert_eq!(dz, 4.0);
    }
    #[test]
    // 少し小さいサイズのモデルで代数的に与えるエッジが矛盾なく繋がるかのテスト
    fn test_linspace_all() {
        let mesh = io::read_vtk("data/Tetra.vtu");
        let ratio = 0.99;
        let cage = CageParameter {
            axis: Vector3::new(0.0, 0.0, 1.0),
            theta0: -PI / 6.0,
            theta1:  PI / 6.0,
            r0: 2.345e-3 * ratio,
            r1: 2.850e-3 * ratio,
            h0: 0.93e-3 * ratio,
            h1: 2.10e-3 * ratio,
            bevel: 0.10e-3 * ratio,
            pocket: PocketParameter {
                x: Vector3::new(0.0, 2.65e-3, 2.00e-3) * ratio,
                r: 0.825e-3 * ratio,
            },
            neck: NeckParameter {
                x: Vector3::new(0.0, 2.65e-3, 2.00e-3 - 0.30e-3) * ratio,
                r: 1.20e-3 * ratio,
                h: 2.45e-3 * ratio,
                dh: 0.1519e-3 * ratio,
                h_ratio: 0.355 * ratio,
                r_ratio: 0.55 * ratio,
            },
        };
        let mut brg = Brg::new(&mesh, &cage);
        let edges = io::read_edge_from_xml("data/face_and_edge_index.xml", "edge").unwrap();
        let faces = io::read_index_from_xml("data/face_and_edge_index.xml", "face").unwrap();    
        brg.set_edge_and_face(&edges, &faces);

        brg.linspace_all(); // ２つ以上の重複定義で，それらが乖離していたらpanicする．
        assert!(true); // panicしなければOK
    }
}




    // pub fn project_all(&mut self) {

    //     let bottom: Vector3<f32>   = Vector3::new(0.0, 0.0, self.cage.h0);
    //     let shoulder: Vector3<f32> = Vector3::new(0.0, 0.0, self.cage.h1);
    //     let top: Vector3<f32>      = Vector3::new(0.0, 0.0, self.cage.neck.h);
    //     let axis: Vector3<f32>     = Vector3::new(0.0, 0.0, 1.0);

    //     for i in self.faces.get("curvature_in").unwrap_or(&HashSet::new()).clone() {
    //         self.mesh.points[i] = self.mesh.points[i].project_on_cylinder(bottom, axis, self.cage.r0);
    //     }
    //     for i in self.faces.get("curvature_out").unwrap_or(&HashSet::new()).clone() {
    //         self.mesh.points[i] = self.mesh.points[i].project_on_cylinder(bottom, axis, self.cage.r1);
    //     }
    //     for i in self.faces.get("bottom").unwrap_or(&HashSet::new()).clone() {
    //         self.mesh.points[i] = self.mesh.points[i].project_on_plane(bottom, axis);
    //     }
    //     for i in self.faces.get("top_left").unwrap_or(&HashSet::new()).clone() {
    //         self.mesh.points[i] = self.mesh.points[i].project_on_plane(top, axis);
    //     }
    //     for i in self.faces.get("top_right").unwrap_or(&HashSet::new()).clone() {
    //         self.mesh.points[i] = self.mesh.points[i].project_on_plane(top, axis);
    //     }
    //     for i in self.faces.get("shoulder_left").unwrap_or(&HashSet::new()).clone() {
    //         self.mesh.points[i] = self.mesh.points[i].project_on_plane(shoulder, axis);
    //     }
    //     for i in self.faces.get("shoulder_right").unwrap_or(&HashSet::new()).clone() {
    //         self.mesh.points[i] = self.mesh.points[i].project_on_plane(shoulder, axis);
    //     }
    //     for i in self.faces.get("sphire").unwrap_or(&HashSet::new()).clone() {
    //         self.mesh.points[i] = self.mesh.points[i].project_on_sphire(self.cage.pocket.x, self.cage.pocket.r);
    //     }
    //     for i in self.faces.get("sphire_left").unwrap_or(&HashSet::new()).clone() {
    //         self.mesh.points[i] = self.mesh.points[i].project_on_sphire(self.cage.neck.x, self.cage.neck.r);
    //     }
    //     for i in self.faces.get("sphire_right").unwrap_or(&HashSet::new()).clone() {
    //         self.mesh.points[i] = self.mesh.points[i].project_on_sphire(self.cage.neck.x, self.cage.neck.r);
    //     }
    // }
