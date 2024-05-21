use std::collections::HashSet;
use std::collections::HashMap;
extern crate nalgebra as na;
use na::Vector3;
use na::Matrix3;

use super::point::Point;
use super::face::Face;
use super::tetra::Tetra;
use super::flower::Flower;

use super::tetra::{tetras_to_face_map, find_surface_faces, find_inner_index, get_neighbor_map, get_surface_map, get_neighbors, make_inverse_map};
use super::math::{get_normal_vector, solve};

fn face_normal(points: &Vec<Point>, face: &[usize; 3]) -> Vector3<f32> {
    let three_points = [points[face[0]].clone(), points[face[1]].clone(), points[face[2]].clone()];
    get_normal_vector(&three_points)
}
fn face_matrix(points: &Vec<Point>, face: &[usize; 3]) -> Matrix3<f32> {
    let mut c: Matrix3<f32> = Matrix3::zeros();
    c.set_column(0, &points[face[1]].clone().as_vec());
    c.set_column(1, &points[face[2]].clone().as_vec());
    c.set_column(2, &points[face[0]].clone().as_vec());
    c
}
fn angular_bisector_matrix(flower: &Flower, points: &Vec<Point>) -> Matrix3<f32> {

    let mut a: Matrix3<f32> = Matrix3::zeros();

    let v0: Vector3<f32> = face_normal(&points, &flower.bottom());
    for i in 0..3 {
        let v1: Vector3<f32> = face_normal(&points, &flower.get(i));
        let bisector: Vector3<f32> = v1 + v0;
        a.set_column(i, &bisector);
    }
    return a
}
pub fn get_intersection_of_flower(flower: &Flower, points: &Vec<Point>) -> Point {
    let a: Matrix3<f32> = angular_bisector_matrix(flower, points);
    let c: Matrix3<f32> = face_matrix(points, &flower.bottom());
    let b: Vector3<f32> = a.component_mul(&c).row_sum_tr();

    Point::from_vec(solve(a.transpose(), b))
}

#[derive(Debug, PartialEq, Clone)]
pub struct Mesh {
    pub points: Vec<Point>,
    pub tetras: Vec<Tetra>,
    pub surface_faces: HashSet<Face>,
    pub inner_index:  HashSet<usize>,
    pub neighbor_map: HashMap<usize, HashSet<usize>>,
    pub surface_map:  HashMap<usize, HashSet<usize>>,
    pub inverse_map:  HashMap<usize, HashSet<Flower>>,
}

impl Mesh {
    pub fn new(points: &Vec<Point>, tetras: &Vec<Tetra>) -> Self {

        let face_map = tetras_to_face_map(&tetras);
        let surface_faces = find_surface_faces(&face_map);
        
        let inner_index: HashSet<usize> = find_inner_index(&face_map, &surface_faces);
        let neighbor_map: HashMap<usize, HashSet<usize>> = get_neighbor_map(&tetras);
        let surface_map: HashMap<usize, HashSet<usize>> = get_surface_map(&surface_faces);
        let inverse_map: HashMap<usize, HashSet<Flower>> = make_inverse_map(&tetras, &inner_index, &face_map);

        Self::load(points, tetras, &surface_faces, &inner_index, &neighbor_map, &surface_map, &inverse_map)
    }

    pub fn save(&self) -> (Vec<Point>, Vec<Tetra>, HashSet<Face>, HashSet<usize>, HashMap<usize, HashSet<usize>>, HashMap<usize, HashSet<usize>>) {
        (self.points.clone(), self.tetras.clone(), self.surface_faces.clone(), self.inner_index.clone(), self.neighbor_map.clone(), self.surface_map.clone())
    }

    pub fn load(points: &Vec<Point>, tetras: &Vec<Tetra>, surface_faces: &HashSet<Face>, inner_index: &HashSet<usize>, neighbor_map: &HashMap<usize, HashSet<usize>>, surface_map: &HashMap<usize, HashSet<usize>>, inverse_map: &HashMap<usize, HashSet<Flower>>) -> Self {
        Self {
            points: points.clone(),
            tetras: tetras.clone(),
            surface_faces: surface_faces.clone(),
            inner_index: inner_index.clone(),
            neighbor_map: neighbor_map.clone(),
            surface_map: surface_map.clone(),
            inverse_map: inverse_map.clone(),
        }
    }
    pub fn smooth_inner(&mut self) {
        let new_points = laplacian_smoothing(&self.points, &self.inner_index, &self.neighbor_map);
        self.points = new_points;
    }
    pub fn smooth_inner_with_cotangent(&mut self) {
        let new_points = cotangent_laplacian_smoothing(&self.points, &self.inner_index, &self.neighbor_map);
        self.points = new_points;
    }
}

pub fn laplacian_smoothing(points: &Vec<Point>, inner_index: &HashSet<usize>, neighbor_map: &HashMap<usize, HashSet<usize>>) -> Vec<Point> {
    let mut new_points = points.clone();

    for &i in inner_index {
        let neighbors = get_neighbors(&neighbor_map, i);
        let sum = neighbors.iter().fold(Point::zero(), |sum, &j| sum.add(new_points[j].as_vec()));
        new_points[i] = sum.mul(1.0 / neighbors.len() as f32);
    }
    new_points
}
pub fn laplacian_smoothing_with_center_normalizing(points: &Vec<Point>, inner_index: &HashSet<usize>, neighbor_map: &HashMap<usize, HashSet<usize>>, center: Vector3<f32>, _: Vector3<f32>, radius: f32) -> Vec<Point> {
    let mut new_points = points.clone();

    for &i in inner_index {
        let neighbors = get_neighbors(&neighbor_map, i);
        let sum = neighbors.iter().fold(Point::zero(), |sum, &j| sum.add(new_points[j].as_vec()));
        let mean = sum.mul(1.0 / neighbors.len() as f32);
        let direction: Vector3<f32> = mean.direction(center);
        new_points[i] = Point::project_on_circle(center, direction, radius);
    }
    new_points
}
pub fn laplacian_smoothing_with_axis_normalizing(points: &Vec<Point>, inner_index: &HashSet<usize>, neighbor_map: &HashMap<usize, HashSet<usize>>, center: Vector3<f32>, axis: Vector3<f32>, radius: f32) -> Vec<Point> {
    let mut new_points = points.clone();

    for &i in inner_index {
        let neighbors = get_neighbors(&neighbor_map, i);
        let sum = neighbors.iter().fold(Point::zero(), |sum, &j| sum.add(new_points[j].as_vec()));
        let mean = sum.mul(1.0 / neighbors.len() as f32);
        let (new_center, direction) = mean.orthogonal(center, axis);
        new_points[i] = Point::project_on_circle(new_center, direction, radius);
    }
    new_points
}
pub fn laplacian_smoothing_with_cone_normalizing(points: &Vec<Point>, inner_index: &HashSet<usize>, neighbor_map: &HashMap<usize, HashSet<usize>>, center: Vector3<f32>, axis: Vector3<f32>, ratio: f32) -> Vec<Point> {

    let mut new_points = points.clone();

    for &i in inner_index {
        let neighbors = get_neighbors(&neighbor_map, i);
        let sum = neighbors.iter().fold(Point::zero(), |sum, &j| sum.add(new_points[j].as_vec()));
        let mean = sum.mul(1.0 / neighbors.len() as f32);
        new_points[i] = mean.orthogonal_cone(center, axis, ratio);
    }
    new_points
}
pub fn laplacian_smoothing_on_plane(points: &Vec<Point>, inner_index: &HashSet<usize>, neighbor_map: &HashMap<usize, HashSet<usize>>, center: Vector3<f32>, normal: Vector3<f32>, _: f32) -> Vec<Point> {
    let mut new_points = points.clone();

    for &i in inner_index {
        let neighbors = get_neighbors(&neighbor_map, i);
        
        let sum = neighbors.iter().fold(Point::zero(), |sum, &j| sum.add(new_points[j].as_vec()));
        let mean = sum.mul(1.0 / neighbors.len() as f32);
        new_points[i] = mean.project_on_plane(center, normal);
    }
    new_points
}
pub fn cotangent_laplacian_smoothing(points: &Vec<Point>, inner_index: &HashSet<usize>, neighbor_map: &HashMap<usize, HashSet<usize>>) -> Vec<Point> {

    let mut new_points = points.clone();

    for &i in inner_index {
        let neighbors = get_neighbors(&neighbor_map, i);
        let mut sum = Vector3::zeros();
        let mut weight = 0.0;
        for &j in neighbors {
            let p0: Vector3<f32> = points[i].as_vec();
            let p1: Vector3<f32> = points[j].as_vec();
            for &k in neighbors {
                if k == j { continue; }
                let p2: Vector3<f32> = points[k].as_vec();
    
                let v0: Vector3<f32> = p0 - p1;
                let v1: Vector3<f32> = p2 - p1;
                let cos_theta = v0.dot(&v1) / (v0.norm() * v1.norm());
                let sin_theta = (1.0 - cos_theta * cos_theta).sqrt().max(1e-3);
                let cot_theta = cos_theta / sin_theta;
                let w = cot_theta / 2.0;
                sum += w * p1;
                weight += w;
            }
        }
        new_points[i] = Point::from_vec(sum / weight);
    }
    new_points
}
pub fn laplacian_smoothing_with_cotangent_and_center_normalizing(points: &Vec<Point>, inner_index: &HashSet<usize>, neighbor_map: &HashMap<usize, HashSet<usize>>, center: Vector3<f32>, _: Vector3<f32>, radius: f32) -> Vec<Point> {
    let mut new_points = points.clone();

    for &i in inner_index {
        let neighbors = get_neighbors(&neighbor_map, i);
        let mut sum: Vector3<f32> = Vector3::zeros();
        let mut weight = 0.0;
        let p0: Vector3<f32> = points[i].as_vec();
        for &j in neighbors {
            let p1: Vector3<f32> = points[j].as_vec();
            for &k in neighbors {
                if k == j { continue; }
                let p2: Vector3<f32> = points[k].as_vec();
    
                let v0: Vector3<f32> = p0 - p1;
                let v1: Vector3<f32> = p2 - p1;
                let cos_theta = v0.dot(&v1) / (v0.norm() * v1.norm());
                let sin_theta = (1.0 - cos_theta * cos_theta).sqrt().max(1e-3);
                let cot_theta = cos_theta / sin_theta;
                let w = cot_theta / 2.0;
                sum += w * p1;
                weight += w;
            }
        }
        let mean = Point::from_vec(sum / weight);
        let direction: Vector3<f32> = mean.direction(center);
        new_points[i] = Point::project_on_circle(center, direction, radius);
    }
    new_points
}
pub fn check_smoothing_quality(old_points: &Vec<Point>, new_points: &Vec<Point>) -> f32 {
    let mut quality = 0.0;
    for (old_point, new_point) in old_points.iter().zip(new_points.iter()) {
        quality += old_point.distance(new_point);
    }
    let n = old_points.len() as f32;
    quality / n
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_face_normal() {
        let points = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
            Point::new(0.0, 1.0, 0.0),
        ];
        let face = [0, 1, 2];
        let normal = face_normal(&points, &face);
        assert_eq!(normal, Vector3::new(0.0, 0.0, 1.0));
    }
    #[test]
    fn test_face_matrix() {
        let points = vec![
            Point::new(1.0,  2.0,  3.0),
            Point::new(10.0, 20.0, 30.0),
            Point::new(100.0, 200.0, 300.0),
        ];
        let face = [0, 1, 2];
        let matrix = face_matrix(&points, &face);
        assert_eq!(matrix, Matrix3::new(
            10.0, 100.0, 1.0,
            20.0, 200.0, 2.0,
            30.0, 300.0, 3.0,
        ));
    }
    #[test]
    fn test_get_intersection_of_flower() {
        let points = vec![
            Point::new( 0.0, 0.0, 3.0),
            Point::new( 3.0, 0.0, 0.0),
            Point::new( 3.0, 3.0, 3.0),
            Point::new( 4.0, 4.0,-1.0),
            Point::new(-1.0, 4.0, 4.0),
            Point::new(-1.0,-1.0,-1.0),
        ];
        let flower = Flower::new(Face::new([0, 1, 2]), [3, 4, 5]);
        let intersection = get_intersection_of_flower(&flower, &points);

        assert!(intersection.distance(&Point::new(0.0, 3.0, 0.0)) < 1e-2);
    }

    #[test]
    fn test_laplacian_smoothing() {
        let sqrt3  = 1.7320508075688772;
        let sqrt27 = sqrt3 * 3.0;

        let points_init = vec![
            Point::new(6.0, 1.0, 0.0),
            Point::new(7.0, 0.0, 1.0),
            Point::new(8.0, -1.0, 0.0),
            Point::new(9.0, 0.0, -1.0),
            Point::new( 6.0, 0.0, 0.0),
            Point::new(-3.0,  sqrt27, 0.0),
            Point::new(-3.0, -sqrt27, 0.0),
        ];
        let inner_index = [0, 1, 2, 3].iter().cloned().collect();
        let neighbor_map = {
            let mut neighbor_map: HashMap<usize, HashSet<usize>> = HashMap::new();
            neighbor_map.insert(0, [1, 2, 3].iter().cloned().collect());
            neighbor_map.insert(1, [0, 4, 5].iter().cloned().collect());
            neighbor_map.insert(2, [0, 4, 6].iter().cloned().collect());
            neighbor_map.insert(3, [0, 5, 6].iter().cloned().collect());
            neighbor_map
        };

        let mut points0: Vec<Point> = points_init.clone();

        let iteration: usize = 50;
        for _ in 0..iteration {
            points0 = laplacian_smoothing(&points0, &inner_index, &neighbor_map);
        }
        assert!(points0[0].distance(&Point::new(0.0, 0.0, 0.0)) < 1e-2);
        assert!(points0[1].distance(&Point::new(1.0, sqrt3, 0.0)) < 1e-2);
        assert!(points0[2].distance(&Point::new(1.0, -sqrt3, 0.0)) < 1e-2);
        assert!(points0[3].distance(&Point::new(-2.0, 0.0, 0.0)) < 1e-2);
    }

    #[test]
    fn test_laplacian_smoothing_with_axis_normalizing() {
        let sqrt3  = 1.7320508075688772;

        let mut points = vec![
            Point::new(0.0, 2.0, 0.0),
            Point::new(0.0, 0.0, 1.0),
            Point::new(0.0, 0.0, 2.0),
            Point::new(2.0, 0.0, 0.0),
        ];
        let inner_index = [1, 2, ].iter().cloned().collect();
        let surface_map = {
            let mut surface_map: HashMap<usize, HashSet<usize>> = HashMap::new();
            surface_map.insert(1, [0, 2].iter().cloned().collect());
            surface_map.insert(2, [1, 3].iter().cloned().collect());
            surface_map
        };
        let iteration: usize = 50;
        for _ in 0..iteration {
            points = laplacian_smoothing_with_axis_normalizing(&points, &inner_index, &surface_map, Vector3::zeros(), Vector3::new(0.0, 0.0, 1.0), 2.0);
        }
        assert_ne!(points[1], Point::new(1.0, sqrt3, 0.0));
        assert_ne!(points[2], Point::new(sqrt3, 1.0, 0.0));
    }
}


