extern crate nalgebra as na;
use std::result::Result;
use std::collections::HashSet;

use na::Vector3;
use na::Matrix3;

use super::point::Point;

pub fn try_solve(a: Matrix3<f32>, b: Vector3<f32>) -> Result<Vector3<f32>, &'static str> {
    let det_a = a.determinant();

    if det_a.abs() < 1e-6 {
        return Err("The matrix is singular and cannot be solved using Cramer's rule.");
    }
    fn get_matrix(a: Matrix3<f32>, b: Vector3<f32>, i: usize) -> f32 {
        let mut a1 = a;
        a1.set_column(i, &b);
        a1.determinant()
    }
    let det_a0 = get_matrix(a, b, 0);
    let det_a1 = get_matrix(a, b, 1);
    let det_a2 = get_matrix(a, b, 2);

    let answer: Vector3<f32> = Vector3::new(det_a0 / det_a, det_a1 / det_a, det_a2 / det_a);
    Ok(answer)
}
pub fn get_normal_vector(points: &[Point; 3]) -> Vector3<f32> {
    let p0:Vector3<f32> = points[0].as_vec();
    let p1:Vector3<f32> = points[1].as_vec();
    let p2:Vector3<f32> = points[2].as_vec();
    let v0:Vector3<f32> = p1 - p0;
    let v1:Vector3<f32> = p2 - p0;
    let normal:Vector3<f32> = v0.cross(&v1);
    let norm = normal.norm();
    normal / norm
}
pub fn sort_three_numbers(v: [usize; 3]) -> ([usize; 3], bool) {
    let a = v[0];
    let b = v[1];
    let c = v[2];

    match (a < b, b < c, a < c) {
        (true, true, true)    => ([a, b, c], true),
        (true, true, false)   => ([a, b, c], true),
        (true, false, true)   => ([a, c, b], false),
        (true, false, false)  => ([c, a, b], true),
        (false, true, true)   => ([b, a, c], false),
        (false, true, false)  => ([b, c, a], true),
        (false, false, true)  => ([b, c, a], true),
        (false, false, false) => ([c, b, a], false),
    }
}
pub fn has_duplication_vec4(index: &[usize; 4]) -> bool {
    let unique_elements: HashSet<_> = index.iter().cloned().collect();
    unique_elements.len() < 4
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_solve() {
        let a = Matrix3::new(
            1.0, 2.0, 3.0,
            4.0, 5.0, 6.0,
            7.0, 8.0, 10.0,
        );
        let b: Vector3<f32> = Vector3::new(1.0, 2.0, 3.0);
        let x: Vector3<f32> = try_solve(a, b).unwrap();
        let y: Vector3<f32> = a.try_inverse().unwrap() * b;
        assert_ne!(x, y);
    }
    #[test]
    fn test_face_sort_three_numbers() {
        assert_eq!(sort_three_numbers([10, 20, 30]), ([10, 20, 30], true));
        assert_eq!(sort_three_numbers([10, 30, 20]), ([10, 20, 30], false));
        assert_eq!(sort_three_numbers([20, 10, 30]), ([10, 20, 30], false));
        assert_eq!(sort_three_numbers([20, 30, 10]), ([10, 20, 30], true));
        assert_eq!(sort_three_numbers([30, 10, 20]), ([10, 20, 30], true));
        assert_eq!(sort_three_numbers([30, 20, 10]), ([10, 20, 30], false));
    }
    #[test]
    fn test_has_duplication_vec4() {
        assert_eq!(has_duplication_vec4(&[1, 2, 3, 4]), false);
        assert_eq!(has_duplication_vec4(&[1, 2, 3, 3]), true);
    }
}
