extern crate nalgebra as na;
use std::result::Result;

use na::Vector3;
use na::Matrix3;

use super::point::Point;

pub fn try_solve(a: Matrix3<f32>, b: Vector3<f32>) -> Result<Vector3<f32>, &'static str> {
    let det_a = a.determinant();

    if det_a.abs() < 1e-2 {
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
pub fn sort_three_numbers(v: [usize; 3]) -> [usize; 3] {
    let a = v[0];
    let b = v[1];
    let c = v[2];
    match (a < b, b < c, a < c) {
        (true, true, _)       => [a, b, c],
        (true, false, true)   => [a, c, b],
        (true, false, false)  => [c, a, b],
        (false, true, true)   => [b, a, c],
        (false, true, false)  => [b, c, a],
        (false, false, true)  => [b, c, a],
        (false, false, false) => [c, b, a],
    }
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
        assert_eq!(sort_three_numbers([1, 2, 3]), [1, 2, 3]);
        assert_eq!(sort_three_numbers([1, 3, 2]), [1, 2, 3]);
        assert_eq!(sort_three_numbers([2, 1, 3]), [1, 2, 3]);
        assert_eq!(sort_three_numbers([2, 3, 1]), [1, 2, 3]);
        assert_eq!(sort_three_numbers([3, 1, 2]), [1, 2, 3]);
        assert_eq!(sort_three_numbers([3, 2, 1]), [1, 2, 3]);
    }
}