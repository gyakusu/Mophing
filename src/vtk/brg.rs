use std::f32::consts::PI;
use std::vec;

use nalgebra as na;
use na::Vector3;

use std::collections::HashMap;
use super::vtk::laplacian_smoothing_on_plane;
use super::vtk::Mesh;
use super::vtk::Face;
use super::vtk::Point;
use super::vtk::laplacian_smoothing_with_axis_normalizing;
use super::vtk::laplacian_smoothing_with_center_normalizing;
// use super::vtk::laplacian_smoothing_with_cone_normalizing;

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
    pub h_ratio: f32, // 面取までの距離／ポケット面までの距離の比率
    pub r_ratio: f32, // 面取半径／ポケット面までの距離の比率
}
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
#[derive(Clone, Debug)]
pub struct Brg {
    mesh: Mesh,
    param: CageParameter,

    edges: HashMap<String, Vec<usize>>,
    faces: HashMap<String, Vec<usize>>,
    periodics: (Vec<usize>, Vec<usize>),
    
}
impl Brg {
    pub fn new(mesh: &Mesh, param: &CageParameter) -> Self {        
        Brg {
            mesh: mesh.clone(),
            param: param.clone(),
            edges: HashMap::new(),   // エッジのインデックス
            faces: HashMap::new(),   // エッジを除く面のインデックス
            periodics: (Vec::new(), Vec::new()), // 周期境界のインデックス
            
        }
    }
    pub fn sample(mesh: &Mesh) -> Self {
        let param = CageParameter::sample();
        Brg::new(mesh, &param)
    }
    pub fn set_edge_and_face(&mut self, edges: HashMap<String, Vec<usize>>, faces: HashMap<String, Vec<usize>>) {
        for (name, index) in edges {
            self.edges.insert(name, index);
        }
        for (name, index) in faces {
            self.faces.insert(name, index);
        }
    }
    pub fn set_periodic(&mut self, left: Vec<usize>, right: Vec<usize>) {
        self.periodics = (left, right);
    }
    pub fn get_points(&self) -> Vec<Point> {
        self.mesh.points.clone()
    }
    pub fn set_points(&mut self, points: Vec<Point>) {
        self.mesh.points = points;
    }
    fn get_curve_params(&self, param: &CageParameter, h: f32, n: usize) -> (Vector3<f32>, f32, f32) {
        let center: Vector3<f32> = Vector3::new(0.0, 0.0, h);
        let theta = param.theta0;
        let dtheta = (param.theta1 - param.theta0) / (n as f32 - 1.0);
        (center, theta, dtheta)
    }
    fn get_asymmetry_curve_params(&self, param: &CageParameter, h: f32, n: usize, r0: f32, r1: f32) -> (Vector3<f32>, f32, f32, f32) {
        let center: Vector3<f32> = Vector3::new(0.0, 0.0, h);
        let y = param.neck.x[1];
        let theta = param.theta0;
        let theta1 = -(calculate_x_diameter(r0, y, r1).0 / r0).asin();
        let dtheta = (theta1 - theta) / (n as f32 - 1.0);
        let radius = r0;
        (center, theta, dtheta, radius)
    }
    fn get_asymmetry_curve_top_params(&self, param: &CageParameter, h: f32, n: usize, r0: f32, r1: f32, r2: f32) -> (Vector3<f32>, f32, f32, f32) {
        let center: Vector3<f32> = Vector3::new(0.0, 0.0, h);
        let y = param.pocket.x[1];
        let theta = (calculate_x_diameter(r0, y, r1).0 / r0).asin();
        let theta1 = (calculate_x_diameter(r0, y, r2).0 / r0).asin();
        let dtheta = (theta1 - theta) / (n as f32 - 1.0);
        let radius = r0;
        (center, theta, dtheta, radius)
    }
    fn get_collar_params(&self, param: &CageParameter, h: f32, n: usize, r: f32, x: f32, y: f32) -> (Vector3<f32>, f32, f32, f32) {
        let center: Vector3<f32> = Vector3::new(x, y, h);
        let y0 = calculate_x_diameter(param.r0, y, r).1;
        let y1 = calculate_x_diameter(param.r1, y, r).1;
        let theta = ((y0 - y) / r).acos();
        let theta1 = ((y1 - y) / r).acos();
        let dtheta = (theta1 - theta) / (n as f32 - 1.0);
        let radius = r;
        (center, -theta, dtheta, radius)
    }
    fn get_ellipse_params(&self, param: &CageParameter, r: f32, n: usize) -> (f32, f32, f32, f32, f32, f32, f32) {
        let x = param.pocket.x[0];
        let y = param.pocket.x[1];
        let z = param.pocket.x[2];
        let dx = calculate_x_diameter(r, y, param.pocket.r).0;
        let dz = calculate_z_diameter(r, y, param.pocket.r);
        let cos_theta = (param.neck.h - param.neck.dh - z) / dz;
        let theta = cos_theta.acos();
        let dtheta = 2.0 * (PI - theta) / (n as f32 - 1.0);        
        (x, z, r, dx, dz, theta, dtheta)
    }
    pub fn get_settings(&self) -> (Vec<Face>, Vec<usize>, HashMap<usize, Vec<usize>>, HashMap<usize, Vec<usize>>) {

        let faces = self.mesh.faces.clone();
        let inner_index = self.mesh.inner_index.clone();
        let neighbor_map = self.mesh.neighbor_map.clone();
        let surface_map = self.mesh.surface_map.clone();

        (faces, inner_index, neighbor_map, surface_map)
    }
    pub fn linspace(&self, name: &str) -> Vec<Point> {
    
        let param = &self.param;
        let n = self.edges.get(name).unwrap().len();

        let r_middle = cutted_sphire_radius(param.pocket.r, param.neck.h - param.pocket.x[2] - param.neck.dh); // 円筒部における半径
        let r_neck   = cutted_sphire_radius(param.neck.r, param.h1 - param.neck.x[2]); // 肩部における爪部の半径
        let r_top    = cutted_sphire_radius(param.neck.r, param.neck.h - param.neck.x[2]); // 頂部における爪外径の半径
        let r_apt    = r_middle + param.neck.dh * param.neck.r_ratio; // 頂部における爪内径の半径

        let mut points: Vec<Point> = Vec::with_capacity(n);

        match name {
            "curve_bottom_in" => {
                let (center, theta, dtheta) = self.get_curve_params(param, param.h0, n);
                (0..n).for_each(|i| points.push(linspace_arc(i, center, param.r0 + param.bevel, theta, dtheta)));
            },
            "curve_bottom_out" => {
                let (center, theta, dtheta) = self.get_curve_params(param, param.h0, n);
                (0..n).for_each(|i| points.push(linspace_arc(i, center, param.r1 - param.bevel, theta, dtheta)));
            },
            "curve_in_bottom" => {
                let (center, theta, dtheta) = self.get_curve_params(param, param.h0 + param.bevel, n);
                (0..n).for_each(|i| points.push(linspace_arc(i, center, param.r0, theta, dtheta)));
            },
            "curve_out_bottom" => {
                let (center, theta, dtheta) = self.get_curve_params(param, param.h0 + param.bevel, n);
                (0..n).for_each(|i| points.push(linspace_arc(i, center, param.r1, theta, dtheta)));
            },
            "curve_in_top_left" | "curve_in_top_right" => {
                let (center, theta, dtheta, radius) = self.get_asymmetry_curve_params(param, param.h1, n, param.r0, r_neck);
                (0..n).for_each(|i| points.push(linspace_arc(i, center, radius, theta, dtheta)));
            },
            "curve_out_top_left" | "curve_out_top_right" => {
                let (center, theta, dtheta, radius) = self.get_asymmetry_curve_params(param, param.h1, n, param.r1, r_neck);
                (0..n).for_each(|i| points.push(linspace_arc(i, center, radius, theta, dtheta)));
            },
            "curve_top_in_left" | "curve_top_in_right" => {
                let (center, theta, dtheta, radius) = self.get_asymmetry_curve_top_params(param, param.neck.h, n, param.r0, r_top, r_apt);
                (0..n).for_each(|i| points.push(linspace_arc(i, center, radius, -theta, -dtheta)));
            }
            "curve_top_out_left" | "curve_top_out_right" => {
                let (center, theta, dtheta, radius) = self.get_asymmetry_curve_top_params(param, param.neck.h, n, param.r1, r_top, r_apt);
                (0..n).for_each(|i| points.push(linspace_arc(i, center, radius, -theta, -dtheta)));
            },
            "straight_left_bottom" | "straight_right_bottom" => {
                let r0 = param.r0 + param.bevel;
                let dr = (param.r1 - param.bevel - r0) / (n as f32 - 1.0);
                let sin = param.theta0.sin();
                let cos = param.theta0.cos();
                let x: Vector3<f32> = Vector3::new(r0 * sin, r0 * cos, param.h0);
                let dx: Vector3<f32> = dr * Vector3::new(sin, cos, 0.0);
                (0..n).for_each(|i| points.push(linspace_line(i, x, dx)));
            },
            "straight_in_left" | "straight_in_right" => {
                let x = param.r0 * param.theta0.sin();
                let y = param.r0 * param.theta0.cos();
                let x0: Vector3<f32> = Vector3::new(x, y, param.h0 + param.bevel);
                let x1: Vector3<f32> = Vector3::new(x, y, param.h1);
                let dx: Vector3<f32> = (x1 - x0) / (n as f32 - 1.0);
                (0..n).for_each(|i| points.push(linspace_line(i, x0, dx)));
            },
            "straight_out_left" | "straight_out_right" => {
                let x = param.r1 * param.theta0.sin();
                let y = param.r1 * param.theta0.cos();
                let x0: Vector3<f32> = Vector3::new(x, y, param.h0 + param.bevel);
                let x1: Vector3<f32> = Vector3::new(x, y, param.h1);
                let dx: Vector3<f32> = (x1 - x0) / (n as f32 - 1.0);
                (0..n).for_each(|i| points.push(linspace_line(i, x0, dx)));
            },
            "straight_top_left" | "straight_top_right" => {
                let sin = param.theta0.sin();
                let cos = param.theta0.cos();
                let (r0, r1) = (param.r0, param.r1);
                let dr = (r1 - r0) / (n as f32 - 1.0);
                let x: Vector3<f32> = Vector3::new(r0 * sin, r0 * cos, param.h1);
                let dx: Vector3<f32> = dr * Vector3::new(sin, cos, 0.0);
                (0..n).for_each(|i| points.push(linspace_line(i, x, dx)));
            },
            "slope_bottom_left_in" | "slope_bottom_right_in" => {
                let sin = param.theta0.sin();
                let cos = param.theta0.cos();
                let (r0, r1) = (param.r0, param.r0 + param.bevel);
                let x0: Vector3<f32> = Vector3::new(r1 * sin, r1 * cos, param.h0);
                let x1: Vector3<f32> = Vector3::new(r0 * sin, r0 * cos, param.h0 + param.bevel);
                let dx: Vector3<f32> = (x1 - x0) / (n as f32 - 1.0);
                (0..n).for_each(|i| points.push(linspace_line(i, x0, dx)));
            }
            "slope_bottom_left_out" | "slope_bottom_right_out" => {
                let sin = param.theta0.sin();
                let cos = param.theta0.cos();
                let(r0, r1) = (param.r1 - param.bevel, param.r1);
                let x0: Vector3<f32> = Vector3::new(r0 * sin, r0 * cos, param.h0);
                let x1: Vector3<f32> = Vector3::new(r1 * sin, r1 * cos, param.h0 + param.bevel);
                let dx: Vector3<f32> = (x1 - x0) / (n as f32 - 1.0);
                (0..n).for_each(|i| points.push(linspace_line(i, x0, dx)));
            }
            "slope_in_top_left" | "slope_in_top_right" => {
                let y = param.neck.x[1];
                let (x0, y0) = calculate_x_diameter(param.r0, y, r_apt);
                let (x1, y1) = calculate_x_diameter(param.r0, y, r_middle);
                let xyz1: Vector3<f32> = Vector3::new(-x0, y0, param.neck.h);
                let xyz0: Vector3<f32> = Vector3::new(-x1, y1, param.neck.h - param.neck.dh * param.neck.h_ratio);
                let dx: Vector3<f32> = (xyz1 - xyz0) / (n as f32 - 1.0);
                (0..n).for_each(|i| points.push(linspace_line(i, xyz0, dx)));
            }
            "slope_out_top_left" | "slope_out_top_right" => {
                let y = param.neck.x[1];
                let (x0, y0) = calculate_x_diameter(param.r1, y, r_apt);
                let (x1, y1) = calculate_x_diameter(param.r1, y, r_middle);
                let xyz1: Vector3<f32> = Vector3::new(-x0, y0, param.neck.h);
                let xyz0: Vector3<f32> = Vector3::new(-x1, y1, param.neck.h - param.neck.dh * param.neck.h_ratio);
                let dx: Vector3<f32> = (xyz1 - xyz0) / (n as f32 - 1.0);
                (0..n).for_each(|i| points.push(linspace_line(i, xyz0, dx)));
            }
            "collar_left_out" | "collar_right_out" => {
                let (center, theta, dtheta, radius) = self.get_collar_params(param, param.h1, n, r_neck, param.neck.x[0], param.neck.x[1]);
                (0..n).for_each(|i| points.push(linspace_arc(i, center, radius, theta, -dtheta)));
            }
            "collar_left_middle" | "collar_right_middle" => {
                let (center, theta, dtheta, radius) = self.get_collar_params(param, param.neck.h, n, r_top, param.neck.x[0], param.neck.x[1]);
                (0..n).for_each(|i| points.push(linspace_arc(i, center, radius, theta, -dtheta)));
            }
            "collar_left_in" | "collar_right_in" => {
                let (center, theta, dtheta, radius) = self.get_collar_params(param, param.neck.h, n, r_apt, param.neck.x[0], param.neck.x[1]);
                (0..n).for_each(|i| points.push(linspace_arc(i, center, radius, theta, -dtheta)));
            },
            "collar_bottom_left" | "collar_bottom_right" => {
                let (center, theta, dtheta, radius) = self.get_collar_params(param, param.neck.h - param.neck.dh, n, r_middle, param.pocket.x[0], param.pocket.x[1]);
                (0..n).for_each(|i| points.push(linspace_arc(i, center, radius, theta, -dtheta)));
            },
            "collar_middle_left" | "collar_middle_right" => {
                let (center, theta, dtheta, radius) = self.get_collar_params(param, param.neck.h - param.neck.dh * param.neck.h_ratio, n, r_middle, param.pocket.x[0], param.pocket.x[1]);
                (0..n).for_each(|i| points.push(linspace_arc(i, center, radius, theta, -dtheta)));
            },
            "ellipse_in_center" => {
                let (x, z, r, dx, dz, theta, dtheta) = self.get_ellipse_params(param, param.r0, n);
                (0..n).for_each(|i| points.push(linspace_ellipse(i, x, z, r, dx, dz, -theta, -dtheta)));
            },
            "ellipse_out_center" => {
                let (x, z, r, dx, dz, theta, dtheta) = self.get_ellipse_params(param, param.r1, n);
                (0..n).for_each(|i| points.push(linspace_ellipse(i, x, z, r, dx, dz, -theta, -dtheta)));
            },
            "ellipse_in_left" | "ellipse_in_right" => {
                let z0 = param.h1;
                let z1 = param.neck.h;
                let dz = (z1 - z0) / (n as f32 - 1.0);
                for i in 0..n {
                    let h = z0 + dz * (i as f32);
                    let r = cutted_sphire_radius(self.param.neck.r, h - param.neck.x[2]);
                    let (x, y) = calculate_x_diameter(self.param.r0, self.param.neck.x[1], r);
                    points.push(Point::new(-x, y, h));
                }
            },
            "ellipse_out_left" | "ellipse_out_right" => {
                let z0 = param.h1;
                let z1 = param.neck.h;
                let dz = (z1 - z0) / (n as f32 - 1.0);
                for i in 0..n {
                    let h = z0 + dz * (i as f32);
                    let r = cutted_sphire_radius(self.param.neck.r, h - param.neck.x[2]);
                    let (x, y) = calculate_x_diameter(self.param.r1, self.param.neck.x[1], r);
                    points.push(Point::new(-x, y, h));
                }
            },
            "stripe_in_left" | "stripe_in_right" => {
                let (x, y) = calculate_x_diameter(param.r0, param.pocket.x[1], r_middle);
                let z0 = param.neck.h - param.neck.dh;
                let z1 = param.neck.h - param.neck.dh * param.neck.h_ratio;
                let dz = (z1 - z0) / (n as f32 - 1.0);
                (0..n).for_each(|i| points.push(Point::new(-x, y, z0 + dz * (i as f32))));
            },
            "stripe_out_left" | "stripe_out_right" => {
                let (x, y) = calculate_x_diameter(param.r1, param.pocket.x[1], r_middle);
                let z0 = param.neck.h - param.neck.dh;
                let z1 = param.neck.h - param.neck.dh * param.neck.h_ratio;
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
    pub fn project_all(&mut self) {

        let bottom: Vector3<f32>   = Vector3::new(0.0, 0.0, self.param.h0);
        let shoulder: Vector3<f32> = Vector3::new(0.0, 0.0, self.param.h1);
        let top: Vector3<f32>      = Vector3::new(0.0, 0.0, self.param.neck.h);
        let axis: Vector3<f32>     = Vector3::new(0.0, 0.0, 1.0);

        for i in self.faces.get("curvature_in").unwrap_or(&Vec::new()).clone() {
            self.mesh.points[i] = self.mesh.points[i].project_on_cylinder(bottom, axis, self.param.r0);
        }
        for i in self.faces.get("curvature_out").unwrap_or(&Vec::new()).clone() {
            self.mesh.points[i] = self.mesh.points[i].project_on_cylinder(bottom, axis, self.param.r1);
        }
        for i in self.faces.get("bottom").unwrap_or(&Vec::new()).clone() {
            self.mesh.points[i] = self.mesh.points[i].project_on_plane(bottom, axis);
        }
        for i in self.faces.get("top_left").unwrap_or(&Vec::new()).clone() {
            self.mesh.points[i] = self.mesh.points[i].project_on_plane(top, axis);
        }
        for i in self.faces.get("top_right").unwrap_or(&Vec::new()).clone() {
            self.mesh.points[i] = self.mesh.points[i].project_on_plane(top, axis);
        }
        for i in self.faces.get("shoulder_left").unwrap_or(&Vec::new()).clone() {
            self.mesh.points[i] = self.mesh.points[i].project_on_plane(shoulder, axis);
        }
        for i in self.faces.get("shoulder_right").unwrap_or(&Vec::new()).clone() {
            self.mesh.points[i] = self.mesh.points[i].project_on_plane(shoulder, axis);
        }
        for i in self.faces.get("sphire").unwrap_or(&Vec::new()).clone() {
            self.mesh.points[i] = self.mesh.points[i].project_on_sphire(self.param.pocket.x, self.param.pocket.r);
        }
        for i in self.faces.get("sphire_left").unwrap_or(&Vec::new()).clone() {
            self.mesh.points[i] = self.mesh.points[i].project_on_sphire(self.param.neck.x, self.param.neck.r);
        }
        for i in self.faces.get("sphire_right").unwrap_or(&Vec::new()).clone() {
            self.mesh.points[i] = self.mesh.points[i].project_on_sphire(self.param.neck.x, self.param.neck.r);
        }
    }
    pub fn smooth_face(&mut self) {

        let points = self.get_points();
        let bottom: Vector3<f32>    = Vector3::new(0.0, 0.0, self.param.h0);
        let shoulder: Vector3<f32>  = Vector3::new(0.0, 0.0, self.param.h1);
        let top: Vector3<f32>       = Vector3::new(0.0, 0.0, self.param.neck.h);
        let axis: Vector3<f32>      = Vector3::new(0.0, 0.0, 1.0);
        let axis_left: Vector3<f32> = Vector3::new(self.param.theta0.cos(), self.param.theta0.sin(), 0.0);
        let axis_right: Vector3<f32> = Vector3::new(self.param.theta1.cos(), self.param.theta1.sin(), 0.0);

        let get_inner_index = |name: &str| self.faces.get(name).unwrap_or(&Vec::new()).clone();
        let surface_map = self.mesh.surface_map.clone();

        let operations: [(&str, fn(&Vec<Point>, Vec<usize>, HashMap<usize, Vec<usize>>, Vector3<f32>, Vector3<f32>, f32) -> Vec<Point>, Vector3<f32>, Vector3<f32>, f32); 12] =
        [
            ("curvature_in", laplacian_smoothing_with_axis_normalizing 
            as fn(&_, _, _, _, _, _) -> _, bottom, axis, self.param.r0),

            ("curvature_out", laplacian_smoothing_with_axis_normalizing 
            as fn(&_, _, _, _, _, _) -> _, bottom, axis, self.param.r1),

            ("bottom", laplacian_smoothing_on_plane 
            as fn(&_, _, _, _, _, _) -> _, bottom, axis, 0.),

            ("top_left", laplacian_smoothing_on_plane 
            as fn(&_, _, _, _, _, _) -> _, top, axis, 0.),

            ("top_right", laplacian_smoothing_on_plane 
            as fn(&_, _, _, _, _, _) -> _, top, axis, 0.),

            ("shoulder_left", laplacian_smoothing_on_plane 
            as fn(&_, _, _, _, _, _) -> _, shoulder, axis, 0.),

            ("shoulder_right", laplacian_smoothing_on_plane 
            as fn(&_, _, _, _, _, _) -> _, shoulder, axis, 0.),

            ("sphire", laplacian_smoothing_with_center_normalizing 
            as fn(&_, _, _, _, _, _) -> _, self.param.pocket.x, axis, self.param.pocket.r),

            ("sphire_left", laplacian_smoothing_with_center_normalizing 
            as fn(&_, _, _, _, _, _) -> _, self.param.neck.x, axis, self.param.neck.r),

            ("sphire_right", laplacian_smoothing_with_center_normalizing 
            as fn(&_, _, _, _, _, _) -> _, self.param.neck.x, axis, self.param.neck.r),

            ("periodic_left", laplacian_smoothing_on_plane 
            as fn(&_, _, _, _, _, _) -> _, bottom, axis_left, 0.),

            ("periodic_right", laplacian_smoothing_on_plane 
            as fn(&_, _, _, _, _, _) -> _, bottom, axis_right, 0.),
        ];
        for (name, function, arg1, arg2, arg3) in operations.iter() {
            let inner_index = get_inner_index(name);
            let smoothed_point = function(&points.clone(), inner_index.clone(), surface_map.clone(), *arg1, *arg2, *arg3);
            inner_index.iter().for_each(|&i| self.mesh.points[i] = smoothed_point[i]);
        }
    }
    pub fn smooth_inner(&mut self) {
        self.mesh.smooth_inner();
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
        let param = CageParameter {
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
        let mut brg = Brg::new(&mesh, &param);
        let edges = io::read_index_from_xml("data/face_and_edge_index.xml", "edge").unwrap();
        let faces = io::read_index_from_xml("data/face_and_edge_index.xml", "face").unwrap();    
        brg.set_edge_and_face(edges, faces);

        brg.linspace_all(); // ２つ以上の重複定義で，それらが乖離していたらpanicする．
        assert!(true); // panicしなければOK
    }
}


