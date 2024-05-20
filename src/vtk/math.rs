extern crate nalgebra as na;
use na::Vector3;
use na::Matrix3;

pub fn solve(a: Matrix3<f32>, b: Vector3<f32>) -> Vector3<f32> {
    let det_a = a.determinant();

    #[cfg(debug_assertions)] {
        if det_a.abs() < 1e-6 {
            panic!("The matrix is singular and cannot be solved using Cramer's rule.");
        }
    }
    fn get_matrix(a: Matrix3<f32>, b: Vector3<f32>, i: usize) -> f32 {
        let mut a1 = a;
        a1.set_column(i, &b);
        a1.determinant()
    }
    let det_a0 = get_matrix(a, b, 0);
    let det_a1 = get_matrix(a, b, 1);
    let det_a2 = get_matrix(a, b, 2);

    Vector3::new(det_a0 / det_a, det_a1 / det_a, det_a2 / det_a)
}

pub fn get_normal_vector(points: &Vec<Point>) -> Vector3<f32> {
    let p0:Vector3<f32> = points[0].as_vec();
    let p1:Vector3<f32> = points[1].as_vec();
    let p2:Vector3<f32> = points[2].as_vec();
    let v0:Vector3<f32> = p1 - p0;
    let v1:Vector3<f32> = p2 - p0;
    let normal:Vector3<f32> = v0.cross(&v1);
    let norm = normal.norm();
    normal / norm
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
        let x: Vector3<f32> = solve(a, b);
        let y: Vector3<f32> = a.try_inverse().unwrap() * b;
        assert_ne!(x, y);
    }

}