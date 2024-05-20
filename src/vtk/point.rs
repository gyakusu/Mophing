extern crate nalgebra as na;
use na::Vector3;

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Point {
    x: Vector3<f32>
}

impl Point {
    pub fn new(x: f32, y: f32, z: f32) -> Self {
        Self {
            x: Vector3::new(x, y, z)
        }
    }
    pub fn from_vec(x: Vector3<f32>) -> Self {
        Self {
            x,
        }
    }
    pub fn zero() -> Self {
        Self {
            x: Vector3::zeros()
        }
    }
    pub fn flip_x(&self) -> Self {
        Self {
            x: Vector3::new(-self.x[0], self.x[1], self.x[2])
        }
    }
    pub fn project_on_circle(center: Vector3<f32>, direction: Vector3<f32>, radius: f32) -> Self {
        #[cfg(debug_assertions)] {
            if (direction.norm_squared() - 1.0).abs() > 1e-6 {
                panic!("direction must be normalized");
            }
        }
        Self {
            x: center + direction * radius
        }
    }
    pub fn project_on_cylinder(&self, center: Vector3<f32>, axis: Vector3<f32>, radius: f32, ) -> Self {
        let (new_center, direction) = self.orthogonal(center, axis);
        Self::project_on_circle(new_center, direction, radius)
    }
    pub fn project_on_sphire(&self, center: Vector3<f32>, radius: f32) -> Self {
        let direction: Vector3<f32> = self.direction(center);
        Self::project_on_circle(center, direction, radius)
    }
    pub fn project_on_plane(&self, center: Vector3<f32>, normal: Vector3<f32>) -> Self {
        #[cfg(debug_assertions)] {
            if (normal.norm_squared() - 1.0).abs() > 1e-6 {
                panic!("normal must be normalized");
            }
        }
        let dx: Vector3<f32> = self.x - center;
        let dot = dx.dot(&normal);
        let new_point: Vector3<f32> = self.x - normal * dot;
        Self {
            x: new_point,
        }
    }
    pub fn as_vec(&self) -> Vector3<f32> {
        self.x
    }
    pub fn add(&self, a: Vector3<f32>) -> Self {
        Self {
            x: self.x + a
        }
    }
    pub fn mul(&self, a: f32) -> Self {
        Self {
            x: self.x * a
        }
    }
    pub fn direction(&self, center: Vector3<f32>) -> Vector3<f32> {
        let dx: Vector3<f32> = self.x - center;
        let distance = dx.norm();
        let direction: Vector3<f32> = if distance != 0.0 {dx / distance} else { Vector3::zeros() };
        direction
    }
    pub fn orthogonal(&self, center: Vector3<f32>, axis: Vector3<f32>) -> (Vector3<f32>, Vector3<f32>) {
        #[cfg(debug_assertions)] {
            if (axis.norm_squared() - 1.0).abs() > 1.0e-6 {
                panic!("axis must be normalized");
            }
        }
        let dx: Vector3<f32> = self.x - center;
        let dot = dx.dot(&axis);
        let new_center: Vector3<f32> = center + axis * dot;
        let orth: Vector3<f32> = dx - axis * dot;
        let distance = orth.norm();
        let direction: Vector3<f32> = if distance != 0.0 {orth / distance} else { Vector3::zeros() };
        (new_center, direction)
    }
    pub fn orthogonal_cone(&self, center: Vector3<f32>, axis: Vector3<f32>, ratio: f32) -> Self {
        #[cfg(debug_assertions)] {
            if (axis.norm_squared() - 1.0).abs() > 1.0e-9 {
                panic!("axis must be normalized");
            }
        }
        let dx: Vector3<f32> = self.x - center;
        let dot = dx.dot(&axis);
        let orth: Vector3<f32> = dx - axis * dot;
        let distance = orth.norm();
        if distance == 0.0 {
            return self.clone();
        }
        let new_length = dot + distance * ratio;
        let new_center: Vector3<f32> = center + axis * new_length;
        let direction: Vector3<f32> = self.x - new_center;
        let new_distance = new_length * ratio / (1. + ratio * ratio);
        let new_direction: Vector3<f32> = direction * new_distance / distance;
        let new_point: Vector3<f32> = new_center + new_direction;
        Self {
            x: new_point,
        }        
    }
    pub fn distance(&self, point: &Self) -> f32 {
        (self.x - point.x).norm()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_direction(){
        let point = Point::new(10.0, 10.0, 10.0);
        let direction: Vector3<f32> = point.direction(Vector3::new(6.0, 7.0, 10.0));
        assert_eq!(direction, Vector3::new(0.8, 0.6, 0.0));
    }
    #[test]
    fn test_orthogonal(){
        let point = Point::new(10.0, 10.0, 10.0);
        let center:Vector3<f32> = Vector3::new(6.0, 7.0, 10.0);
        let axis:Vector3<f32> = Vector3::new(1.0, 0.0, 0.0);
        let (new_center, direction) = point.orthogonal(center, axis);
        assert_eq!(new_center, Vector3::new(10.0, 7.0, 10.0));
        assert_eq!(direction, Vector3::new(0.0, 1.0, 0.0));
    }
    #[test]
    fn test_orthogonal_cone(){
        let sqrt1_3  = 1.7320508075688772 / 3.0;
        let point = Point::new(10.0, 10.0, 10.0);
        let center0:Vector3<f32> = Vector3::new(9.0, 10.0 - sqrt1_3, 10.0);
        let center1:Vector3<f32> = Vector3::new(10.0, 10.0, 10.0);
        let center2:Vector3<f32> = Vector3::new(5.0, 10.0, 10.0);
        let center3:Vector3<f32> = Vector3::new(9.0, 10.0, 7.0);
        let axis:Vector3<f32> = Vector3::new(1.0, 0.0, 0.0);

        let new_point0 = point.orthogonal_cone(center0, axis, sqrt1_3);
        assert_eq!(new_point0, point);

        let new_point1 = point.orthogonal_cone(center1, axis, sqrt1_3);
        assert_eq!(new_point1, point);

        let new_point2 = point.orthogonal_cone(center2, axis, sqrt1_3);
        assert_eq!(new_point2, point);

        let new_point3 = point.orthogonal_cone(center3, axis, 1.0);
        assert_eq!(new_point3, Point::new(11.0, 10.0, 9.0));

    }
    #[test]
    fn test_project_on_cylinder() {
        let point = Point::new(-8.0, -6.0, 1.0);
        let center:Vector3<f32> = Vector3::zeros();
        let axis:Vector3<f32> = Vector3::new(0.0, 0.0, 1.0);
        let radius = 5.0;
        let new_point = point.project_on_cylinder(center, axis, radius);
        assert_eq!(new_point, Point::new(-4.0, -3.0, 1.0));
    }
    #[test]
    fn test_project_on_sphire() {
        let point = Point::new(-8.0, -6.0, 0.0);
        let center:Vector3<f32> = Vector3::zeros();
        let radius = 5.0;
        let new_point = point.project_on_sphire(center, radius);
        assert_eq!(new_point, Point::new(-4.0, -3.0, 0.0));
    }

}
