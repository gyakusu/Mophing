use core::panic;
use std::collections::HashSet;
use std::collections::HashMap;
use std::result::Result;
use std::vec;
extern crate nalgebra as na;
use na::Vector3;
use na::Matrix3;

use super::point::Point;
use super::face::Face;
use super::tetra::Tetra;
use super::flower::Flower;

use super::tetra::{tetras_to_face_map, find_surface_faces, find_inner_index, get_neighbor_map, get_surface_map, get_neighbors, make_inverse_map};
use super::math::{get_normal_vector, try_solve};

fn face_normal(points: &Vec<Point>, face: &[usize; 3]) -> Vector3<f32> {
    let three_points = [points[face[0]].clone(), points[face[1]].clone(), points[face[2]].clone()];
    get_normal_vector(&three_points)
}
fn face_matrix(points: &Vec<Point>, face: &[usize; 3]) -> Matrix3<f32> {
    let mut c: Matrix3<f32> = Matrix3::zeros();
    c.set_row(0, &points[face[1]].as_vec().transpose());
    c.set_row(1, &points[face[2]].as_vec().transpose());
    c.set_row(2, &points[face[0]].as_vec().transpose());
    c
}
fn tetra_volume(points: &Vec<Point>, tetra: &Tetra) -> f32 {
    let t = tetra.index;
    let a: Vector3<f32> = points[t[1]].as_vec() - points[t[0]].as_vec();
    let b: Vector3<f32> = points[t[2]].as_vec() - points[t[0]].as_vec();
    let c: Vector3<f32> = points[t[3]].as_vec() - points[t[0]].as_vec();
    let volume = a.dot(&b.cross(&c)) / 6.0;
    volume
}
fn tetra_area(points: &Vec<Point>, tetra: &Tetra) -> Vec<f32> {
    let t = tetra.index;
    let a: Vector3<f32> = points[t[0]].as_vec();
    let b: Vector3<f32> = points[t[1]].as_vec();
    let c: Vector3<f32> = points[t[2]].as_vec();
    let d: Vector3<f32> = points[t[3]].as_vec();

    let area = vec![
        (d - b).cross(&(c - b)).norm() / 2.0,
        (d - c).cross(&(a - c)).norm() / 2.0,
        (d - a).cross(&(b - a)).norm() / 2.0,
        (c - b).cross(&(a - b)).norm() / 2.0,
    ];
    area
}
fn sort_vec_index(v: &Vec<f32>) -> Vec<usize> {
    let mut array: [usize; 4] = [0, 1, 2, 3];
    array.sort_by(|&a, &b| v[a].partial_cmp(&v[b]).unwrap());
    array.to_vec()
}
fn unit_normal(points: &Vec<Point>, tetra: &Tetra, i: usize) -> Vector3<f32> {
    let t = tetra.index;
    let a: Vector3<f32> = points[t[1]].as_vec() - points[t[0]].as_vec();
    let b: Vector3<f32> = points[t[2]].as_vec() - points[t[1]].as_vec();
    let c: Vector3<f32> = points[t[3]].as_vec() - points[t[2]].as_vec();
    let d: Vector3<f32> = points[t[0]].as_vec() - points[t[3]].as_vec();
    let normal = match i {
        0 => c.cross(&b),
        1 => c.cross(&d),
        2 => a.cross(&d),
        3 => a.cross(&b),
        _ => Vector3::zeros(),
    };
    normal.normalize()
}
fn angular_bisector_matrix(flower: &Flower, points: &Vec<Point>) -> Matrix3<f32> {
    let v: Vector3<f32> = face_normal(&points, &flower.bottom());

    let mut a: Matrix3<f32> = Matrix3::zeros();
    for i in 0..3 {
        let f0 = flower.get(i);
        let v1: Vector3<f32> = face_normal(&points, &f0);

        let bisector: Vector3<f32> = v1 + v;
        let bisector_norm = bisector.norm();

        let b: Vector3<f32> =  if bisector_norm < 1e-3 {
            let v2: Vector3<f32> = points[f0[1]].as_vec() - points[f0[0]].as_vec();
            let v3: Vector3<f32> = v.cross(&v2).normalize();
            v3
        }
        else {
            let distance = v.dot(&(points[f0[0]].as_vec() - points[f0[2]].as_vec()));
            let distance_fliped = distance * flower.flip();
    
            let is_flip = if distance_fliped > 0.0 {1.0} else {-1.0};
            bisector / bisector_norm * is_flip
        };
        a.set_row(i, &b.transpose());     
    }
    a
}
fn angular_quadrisector_matrix(flower: &Flower, points: &Vec<Point>) -> Matrix3<f32> {
    let v: Vector3<f32> = face_normal(&points, &flower.bottom());

    let mut a: Matrix3<f32> = Matrix3::zeros();
    for i in 0..3 {
        let f0 = flower.get(i);
        let v1: Vector3<f32> = face_normal(&points, &f0);

        let distance = v.dot(&(points[f0[0]].as_vec() - points[f0[2]].as_vec()));
        let distance_fliped = distance * flower.flip();
        let is_flip = if distance_fliped > 0.0 {1.0} else {-1.0};
        let v2: Vector3<f32> = points[f0[1]].as_vec() - points[f0[0]].as_vec();
        let v3: Vector3<f32> = v.cross(&v2).normalize();

        let bisector: Vector3<f32> = v1 + v;
        let bisector_norm = bisector.norm();

        let b: Vector3<f32> =  if bisector_norm < 1e-3 {
            v3
        }
        else {
            bisector / bisector_norm * is_flip
        }; 
        let bisector: Vector3<f32> = b + v;
        let bisector_norm = bisector.norm();

        let b: Vector3<f32> =  if bisector_norm < 1e-3 {
            v3
        }
        else {
            bisector / bisector_norm * is_flip
        };

        let bisector: Vector3<f32> = b + v;
        let bisector_norm = bisector.norm();

        let b: Vector3<f32> =  if bisector_norm < 1e-3 {
            v3
        }
        else {
            bisector / bisector_norm * is_flip
        };
        a.set_row(i, &b.transpose());
    }
    a
}
pub fn intersect_flower(flower: &Flower, points: &Vec<Point>, use_quad: bool) -> Result<Point, &'static str> {

    let a: Matrix3<f32> = if use_quad { 
        angular_quadrisector_matrix(flower, points)
    }
    else {
        angular_bisector_matrix(flower, points)
    };
    let c: Matrix3<f32> = face_matrix(points, &flower.bottom());
    let b: Vector3<f32> = a.component_mul(&c).column_sum();

    let d: Result<Vector3<f32>, &'static str> = try_solve(a, b);
    match d {
        Ok(answer) => {
            Ok(Point::from_vec(answer))
        },
        Err(_) => Err("The matrix is singular and cannot be solved using Cramer's rule."),
    }
}
pub fn flower_smoothing(points: &Vec<Point>, inner_index: &HashSet<usize>, inverse_map: &HashMap<usize, HashSet<Flower>>) -> Vec<Point> {
    let mut new_points = points.clone();

    for &i in inner_index {
        let flowers = inverse_map.get(&i).unwrap();
        let mut sum = Vector3::zeros();
        let mut waight = 0.0;
        for flower in flowers {
            let intersection = intersect_flower(&flower, &points, true);
            match intersection {
                Ok(p) => {
                    sum += p.as_vec();
                    waight += 1.0;
                },
                Err(_) => {
                },
            }
        }
        if waight < 1e-1 {
            continue;
        }
        new_points[i] = Point::from_vec(sum / waight);
    }
    new_points
}
pub fn flip_negative_volume(points: &Vec<Point>, tetras: &Vec<Tetra>, inner_index: &HashSet<usize>, ratio: f32, use_sort: bool) -> (Vec<Point>, bool) {

    let mut new_points = points.clone();
    let mut num_fliped = 0;
    for tetra in tetras {
        let volume = tetra_volume(&new_points, &tetra);
        if volume > 0.0 {
            continue;
        }
        num_fliped += 1;
        let t = tetra.index;
        let area = tetra_area(&new_points, &tetra);
        // let index = sort_vec_index(&area);
        let index = if use_sort { sort_vec_index(&area) } else { [0, 1, 2, 3].to_vec() };
        let i = if inner_index.contains(&t[index[0]]) { index[0] }
        else if inner_index.contains(&t[index[1]]) { index[1] }
        else if inner_index.contains(&t[index[2]]) { index[2] }
        else if inner_index.contains(&t[index[3]]) { index[3] }
        else {
            println!("tetras num {}", tetras.len());
            println!("first tetra {:?}", tetras[0].index);
            println!("now tetra {:?}", t);
            println!(". The tetra on {}, {}, {}, {} is error. ", t[index[0]], t[index[1]], t[index[2]], t[index[3]]);
            panic!("The tetra is on the surface! Can't flip it!");
        };
        print!("{:?}", i);
        let height = volume * 3.0 / area[i]; // < 0
        let length = f32::max(height * (-1.0 - ratio), 1e-6);
        let normal: Vector3<f32> = unit_normal(&new_points, &tetra, i);
        new_points[t[i]] = new_points[t[i]].add(normal * length);
    }
    let success = num_fliped == 0;
    if success {
        print!("There are no negative volume tetras.");
    }
    println!("");
    (new_points, success)
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
pub fn normalize_center(points: &Vec<Point>, center: Vector3<f32>, radius: f32, inner_index: &HashSet<usize>) -> Vec<Point> {
    let mut new_points = points.clone();

    for &i in inner_index {
        let direction: Vector3<f32> = new_points[i].direction(center);
        new_points[i] = Point::project_on_circle(center, direction, radius);
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
pub fn laplacian_smoothing_wrapper(points: &Vec<Point>, inner_index: &HashSet<usize>, neighbor_map: &HashMap<usize, HashSet<usize>>, _: Vector3<f32>, _: Vector3<f32>, _: f32) -> Vec<Point> {
    laplacian_smoothing(points, inner_index, neighbor_map)
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
pub fn reindex_surface(surface_faces: &HashSet<(Face, bool)>) -> HashMap<usize, usize> {
    let mut index_map: HashMap<usize, usize> = HashMap::new();
    let mut new_index = 0;

    for face in surface_faces {
        for i in 0..3 {
            let index = face.0.get(i);
            if !index_map.contains_key(&index) {
                index_map.insert(index, new_index);
                new_index += 1;
            }
        }
    }
    index_map
}

#[derive(Debug, PartialEq, Clone)]
pub struct Mesh {
    pub points: Vec<Point>,
    pub tetras: Vec<Tetra>,
    pub surface_faces: HashSet<(Face, bool)>,
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
    pub fn save(&self) -> (Vec<Point>, Vec<Tetra>, HashSet<(Face, bool)>, HashSet<usize>, HashMap<usize, HashSet<usize>>, HashMap<usize, HashSet<usize>>) {
        (self.points.clone(), self.tetras.clone(), self.surface_faces.clone(), self.inner_index.clone(), self.neighbor_map.clone(), self.surface_map.clone())
    }
    pub fn load(points: &Vec<Point>, tetras: &Vec<Tetra>, surface_faces: &HashSet<(Face, bool)>, inner_index: &HashSet<usize>, neighbor_map: &HashMap<usize, HashSet<usize>>, surface_map: &HashMap<usize, HashSet<usize>>, inverse_map: &HashMap<usize, HashSet<Flower>>) -> Self {
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
    pub fn smooth_inner_with_flower(&mut self) {
        let new_points = flower_smoothing(&self.points, &self.inner_index, &self.inverse_map);
        self.points = new_points;
    }
    pub fn flip_negative_volume(&mut self, inner_index: &HashSet<usize>, ratio: f32, use_sort: bool) -> bool {
        let (new_points, success) = flip_negative_volume(&self.points, &self.tetras, inner_index, ratio, use_sort);
        self.points = new_points;
        success
    }
    pub fn normalize_center(&mut self, center: Vector3<f32>, radius: f32, inner_index: &HashSet<usize>) {
        let new_points = normalize_center(&self.points, center, radius, inner_index);
        self.points = new_points;
    }
    pub fn extract_surface(&self) -> (HashMap<usize, usize>, HashSet<(Face, bool)>) {
        (reindex_surface(&self.surface_faces), self.surface_faces.clone())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const TETRAS0 : [Tetra; 8] = [
        Tetra { index: [0, 1, 3, 5] },
        Tetra { index: [0, 3, 2, 5] },
        Tetra { index: [0, 2, 4, 5] },
        Tetra { index: [0, 4, 1, 5] },
        Tetra { index: [0, 1, 4, 6] },
        Tetra { index: [0, 4, 2, 6] },
        Tetra { index: [0, 2, 3, 6] },
        Tetra { index: [0, 3, 1, 6] },
    ];
    const POINTS0 : [Vector3<f32>; 7] = [
        Vector3::new( 0.0, 0.0, 0.0),
        Vector3::new( 1.0, 0.0, 0.0),
        Vector3::new( -1.0, 0.0, 0.0),
        Vector3::new( 0.0, 1.0, 0.0),
        Vector3::new( 0.0, -1.0, 0.0),
        Vector3::new( 0.0, 0.0, 1.0),
        Vector3::new( 0.0, 0.0, -1.0),
    ];
    const INNER_INDEX0: [usize; 1] = [0, ];

    #[test]
    fn test_reindex_surface() {
        let tetras: Vec<Tetra> = TETRAS0.to_vec();
        let face_map = tetras_to_face_map(&tetras);
        let surface_faces = find_surface_faces(&face_map);
        let index_map = reindex_surface(&surface_faces);
        assert!(index_map.len() == 6);
        assert!(index_map.contains_key(&1));
        assert!(index_map.contains_key(&6));
        assert!(!index_map.contains_key(&0));
    }
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
            Point::new(100.0, 200.0, 300.0),
            Point::new(1.0,  2.0,  3.0),
            Point::new(10.0, 20.0, 30.0),
        ];
        let face = [0, 1, 2];
        let matrix = face_matrix(&points, &face);
        assert_eq!(matrix, Matrix3::new(
            1.0,  2.0,  3.0,
            10.0, 20.0, 30.0,
            100.0, 200.0, 300.0,
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
        let flower = Flower::from_vec([0, 1, 2, 3, 4, 5]);
        let intersection = intersect_flower(&flower, &points, false).unwrap();
        let p = intersection;

        assert!(p.distance(&Point::new(0.0, 3.0, 0.0)) < 1e-2);
    }
    #[test]
    fn test_flower_smoothing() {
        let points: Vec<Point> = (0..7).map(|i| Point::from_vec(POINTS0[i])).collect();
        let tetras: Vec<Tetra> = TETRAS0.to_vec();
        let inner_index: HashSet<usize> = INNER_INDEX0.iter().cloned().collect();
        let face_map = tetras_to_face_map(&tetras);
        let inverse_map = make_inverse_map(&tetras, &inner_index, &face_map);

        let mut points0: Vec<Point> = points.clone();
        points0 = flower_smoothing(&points0, &inner_index, &inverse_map);
        assert!(points0[0].distance(&Point::new(0.0, 0.0, 0.0)) < 1e-2);
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
        // assert_ne!(points[1], Point::new(1.0, sqrt3, 0.0));
        // assert_ne!(points[2], Point::new(sqrt3, 1.0, 0.0));
        assert!(points[1].distance(&Point::new(1.0, sqrt3, 0.0)) < 1e-2);
        assert!(points[2].distance(&Point::new(sqrt3, 1.0, 0.0)) < 1e-2);
    }
    #[test]
    fn test_volume() {
        let points = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
            Point::new(0.0, 1.0, 0.0),
            Point::new(0.0, 0.0, 1.0),
        ];
        let tetra = Tetra { index: [0, 1, 2, 3] };
        let volume = tetra_volume(&points, &tetra);
        assert!(volume - 1.0 / 6.0 < 1e-6);
    }
    #[test]
    fn test_area() {
        let sqrt3 = 1.7320508075688772;
        let points = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
            Point::new(0.0, 1.0, 0.0),
            Point::new(0.0, 0.0, 1.0),
        ];
        let tetra = Tetra { index: [0, 1, 2, 3] };
        let area = tetra_area(&points, &tetra);
        assert_eq!(area[0], sqrt3 / 2.0);
        assert_eq!(area[1], 0.5);
        assert_eq!(area[2], 0.5);
        assert_eq!(area[3], 0.5);
    }
    #[test]
    fn test_sort_vec_index() {
        let v = vec![3.0, 1.0, 2.0, 0.0];
        let index = sort_vec_index(&v);
        assert_eq!(index, vec![3, 1, 2, 0]);
    }
}


