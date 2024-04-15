use std::collections::HashSet;
use std::collections::HashMap;
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
            if (direction.norm() - 1.0).abs() > 1e-6 {
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
            if (normal.norm() - 1.0).abs() > 1e-6 {
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
            if (axis.norm() - 1.0).abs() > 1.0e-6 {
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
            if (axis.norm() - 1.0).abs() > 1.0e-9 {
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
}

#[derive(Debug, Eq, Hash, PartialEq, Clone, Copy)]
pub struct Face {
    pub index: [usize; 3],
}

impl Face {
    pub fn new(index: [usize; 3]) -> Result<Self, &'static str> {
        if index.len() != 3 {
            return Err("index must have a length of 3");
        }        
        let mut sorted_index = index;
        sorted_index.sort();
        Ok(Self {
            index: sorted_index,
        })
    }
    pub fn as_i64(&self) -> [usize; 3] {
        self.index
    }
}

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Tetra {
    pub index: [usize; 4],
}
impl Tetra {
    pub fn new(index: [usize; 4]) -> Result<Self, &'static str> {
        if index.len() != 4 {
            return Err("index must have a length of 4");
        }        
        let mut sorted_index = index;
        sorted_index.sort();
        Ok(Self {
            index: sorted_index,
        })
    }
}

#[derive(Debug, PartialEq, Clone)]
pub struct Mesh {
    pub points: Vec<Point>,
    pub tetras: Vec<Tetra>,
    pub faces: Vec<Face>,
    pub outer_index: Vec<usize>,
    pub inner_index: Vec<usize>,
    pub neighbor_map: HashMap<usize, Vec<usize>>,
    pub outer_map: HashMap<usize, Vec<usize>>,
}

impl Mesh {
    pub fn new(points: Vec<Point>, tetras: Vec<Tetra>) -> Self {
        
        let faces: Vec<Face> = tetras_to_faces(&tetras);
        let outer_index: Vec<usize> = find_outer_index(faces.clone());
        let inner_index: Vec<usize> = find_inner_index(faces.clone(), outer_index.clone());
        let neighbor_map: HashMap<usize, Vec<usize>> = find_neighbors(&tetras);
        let outer_map: HashMap<usize, Vec<usize>> = find_outer_neighbors(&neighbor_map, &inner_index);

        Self::load(points, tetras, faces, outer_index, inner_index, neighbor_map, outer_map)
    }

    pub fn save(&self) -> (Vec<Point>, Vec<Tetra>, Vec<Face>, Vec<usize>, Vec<usize>, HashMap<usize, Vec<usize>>, HashMap<usize, Vec<usize>>) {
        (self.points.clone(), self.tetras.clone(), self.faces.clone(), self.outer_index.clone(), self.inner_index.clone(), self.neighbor_map.clone(), self.outer_map.clone())
    }

    pub fn load(points: Vec<Point>, tetras: Vec<Tetra>, faces: Vec<Face>, outer_index: Vec<usize>, inner_index: Vec<usize>, neighbor_map: HashMap<usize, Vec<usize>>, outer_map: HashMap<usize, Vec<usize>>) -> Self {
        Self {
            points,
            tetras,
            faces,
            outer_index,
            inner_index,
            neighbor_map,
            outer_map,
        }
    }
    pub fn smooth_inner(&mut self) {
        let new_points = laplacian_smoothing(self.points.clone(), self.inner_index.clone(), self.neighbor_map.clone());
        self.points = new_points;
    }
}
pub fn tetras_to_faces(tetras: &Vec<Tetra>) -> Vec<Face> {
    let mut faces: Vec<Face> = Vec::new();
    for tetra in tetras {
        for i in 0..4 {
            for j in (i + 1)..4 {
                for k in (j + 1)..4 {
                    let face = Face {
                        index: [tetra.index[i], tetra.index[j], tetra.index[k]],
                    };
                    faces.push(face);
                }
            }
        }
    }
    faces
}
pub fn find_outer_index(faces: Vec<Face>) -> Vec<usize> {

    let mut face_counts: HashMap<Face, usize> = HashMap::new();
    for face in &faces {
        *face_counts.entry(face.clone()).or_insert(0) += 1;
    }
    let shared_faces: HashSet<Face> = face_counts.iter()
        .filter(|&(_face, &count)| count > 1)
        // .filter(|&(_face, &count)| count == 2)
        .map(|(face, _count)| face.clone())
        .collect();
    let not_shared_faces: Vec<Face> = faces.iter()
        .filter(|face| !shared_faces.contains(face))
        .cloned()
        .collect();
    let outer_map: HashSet<usize> = not_shared_faces.iter()
        .flat_map(|face| face.index.iter())
        .cloned()
        .collect();
    let mut outer_index: Vec<usize> = outer_map.iter().cloned().collect();
    outer_index.sort();
    outer_index
}
pub fn find_inner_index(faces: Vec<Face>, outer_index: Vec<usize>) -> Vec<usize> {

    let all_index: HashSet<usize> = faces.iter()
        .flat_map(|face| face.index.iter())
        .cloned()
        .collect();
    let mut inner_index: Vec<usize> = all_index.difference(&outer_index.iter().cloned().collect())
        .cloned()
        .collect();
    inner_index.sort();
    inner_index
}
pub fn find_neighbors(tetras: &Vec<Tetra>) -> HashMap<usize, Vec<usize>> {
    let mut neighbor_map: HashMap<usize, Vec<usize>> = HashMap::new();
    for tetra in tetras {
        for i in 0..4 {
            for j in 0..4 {
                if i != j {
                    neighbor_map.entry(tetra.index[i]).or_insert_with(Vec::new).push(tetra.index[j]);
                }
            }
        }
    }
    for (_key, value) in neighbor_map.iter_mut() {
        value.sort_unstable(); // ソート
        value.dedup(); // 重複の削除
    }
    neighbor_map
}
pub fn find_outer_neighbors(neighbor_map: &HashMap<usize, Vec<usize>>, inner_index: &[usize]) -> HashMap<usize, Vec<usize>> {
    let mut outer_map: HashMap<usize, Vec<usize>> = HashMap::new();
    for (key, value) in neighbor_map {
        let filtered_value: Vec<usize> = value.iter().filter(|&index| !inner_index.contains(index)).cloned().collect();
        outer_map.insert(*key, filtered_value);
    }
    outer_map
}

pub fn laplacian_smoothing(points: Vec<Point>, inner_index: Vec<usize>, neighbor_map: HashMap<usize, Vec<usize>>) -> Vec<Point> {
    let mut new_points = points.clone();

    for &i in &inner_index {
        let neighbors = &neighbor_map[&i];
        let sum = neighbors.iter().fold(Point::zero(), |sum, &j| sum.add(new_points[j].as_vec()));
        new_points[i] = sum.mul(1.0 / neighbors.len() as f32);
    }
    new_points
}
pub fn laplacian_smoothing_with_center_normalizing(points: Vec<Point>, inner_index: Vec<usize>, neighbor_map: HashMap<usize, Vec<usize>>, center: Vector3<f32>, _: Vector3<f32>, radius: f32) -> Vec<Point> {
    let mut new_points = points.clone();

    for &i in &inner_index {
        let neighbors = &neighbor_map[&i];
        let sum = neighbors.iter().fold(Point::zero(), |sum, &j| sum.add(new_points[j].as_vec()));
        let mean = sum.mul(1.0 / neighbors.len() as f32);
        let direction: Vector3<f32> = mean.direction(center);
        new_points[i] = Point::project_on_circle(center, direction, radius);
    }
    new_points
}
pub fn laplacian_smoothing_with_axis_normalizing(points: Vec<Point>, inner_index: Vec<usize>, neighbor_map: HashMap<usize, Vec<usize>>, center: Vector3<f32>, axis: Vector3<f32>, radius: f32) -> Vec<Point> {
    let mut new_points = points.clone();

    for &i in &inner_index {
        let neighbors = &neighbor_map[&i];
        let sum = neighbors.iter().fold(Point::zero(), |sum, &j| sum.add(new_points[j].as_vec()));
        let mean = sum.mul(1.0 / neighbors.len() as f32);
        let (new_center, direction) = mean.orthogonal(center, axis);
        new_points[i] = Point::project_on_circle(new_center, direction, radius);
    }
    new_points
}
pub fn laplacian_smoothing_with_cone_normalizing(points: Vec<Point>, inner_index: Vec<usize>, neighbor_map: HashMap<usize, Vec<usize>>, center: Vector3<f32>, axis: Vector3<f32>, ratio: f32) -> Vec<Point> {

    let mut new_points = points.clone();

    for &i in &inner_index {
        let neighbors = &neighbor_map[&i];
        let sum = neighbors.iter().fold(Point::zero(), |sum, &j| sum.add(new_points[j].as_vec()));
        let mean = sum.mul(1.0 / neighbors.len() as f32);
        new_points[i] = mean.orthogonal_cone(center, axis, ratio);
    }
    new_points
}
pub fn laplacian_smoothing_on_plane(points: Vec<Point>, inner_index: Vec<usize>, neighbor_map: HashMap<usize, Vec<usize>>, center: Vector3<f32>, normal: Vector3<f32>, _: f32) -> Vec<Point> {
    let mut new_points = points.clone();

    for &i in &inner_index {
        let neighbors = &neighbor_map[&i];
        let sum = neighbors.iter().fold(Point::zero(), |sum, &j| sum.add(new_points[j].as_vec()));
        let mean = sum.mul(1.0 / neighbors.len() as f32);
        new_points[i] = mean.project_on_plane(center, normal);
    }
    new_points
}
#[cfg(test)]
mod tests {
    use super::*;

    const TETRAS0 : [Tetra; 4] = [
        Tetra { index: [0, 1, 2, 3] },
        Tetra { index: [0, 1, 2, 4] },
        Tetra { index: [0, 1, 3, 4] },
        Tetra { index: [0, 2, 3, 4] },
    ];
    const TETRAS1 : [Tetra; 8] = [
        Tetra { index: [0, 1, 3, 5] },
        Tetra { index: [0, 1, 3, 6] },
        Tetra { index: [0, 1, 4, 5] },
        Tetra { index: [0, 1, 4, 6] },
        Tetra { index: [0, 2, 3, 5] },
        Tetra { index: [0, 2, 3, 6] },
        Tetra { index: [0, 2, 4, 5] },
        Tetra { index: [0, 2, 4, 6] },
    ];

    #[test]
    fn test_find_outer_index() {
        let tetras = TETRAS0.to_vec();
        let faces = tetras_to_faces(&tetras);
        let outer_index = find_outer_index(faces);
        assert_eq!(outer_index, vec![1, 2, 3, 4]);
    }
    #[test]
    fn test_find_inner_index() {
        let tetras = TETRAS0.to_vec();
        let faces = tetras_to_faces(&tetras);
        let outer_index = vec![1, 2, 3, 4];
        let find_inner_index = find_inner_index(faces.clone(), outer_index);
        assert_eq!(find_inner_index, vec![0, ]);
    }
    #[test]
    fn test_find_neighbors() {
        let tetras = TETRAS1.to_vec();
        let neighbor_map = find_neighbors(&tetras);
        assert_eq!(neighbor_map.get(&0), Some(&vec![1, 2, 3, 4, 5, 6]));
        assert_eq!(neighbor_map.get(&1), Some(&vec![0, 3, 4, 5, 6]));
        assert_eq!(neighbor_map.get(&2), Some(&vec![0, 3, 4, 5, 6]));
        assert_eq!(neighbor_map.get(&3), Some(&vec![0, 1, 2, 5, 6]));
        assert_eq!(neighbor_map.get(&4), Some(&vec![0, 1, 2, 5, 6]));
        assert_eq!(neighbor_map.get(&5), Some(&vec![0, 1, 2, 3, 4]));
        assert_eq!(neighbor_map.get(&6), Some(&vec![0, 1, 2, 3, 4]));
    }
    #[test]
    fn test_find_outer_neighbors() {
        let tetras = TETRAS1.to_vec();
        let neighbor_map = find_neighbors(&tetras);
        let inner_index = vec![0, 6];
        let outer_map = find_outer_neighbors(&neighbor_map, &inner_index);
        assert_eq!(outer_map.get(&1), Some(&vec![3, 4, 5, ]));
        assert_eq!(outer_map.get(&2), Some(&vec![3, 4, 5, ]));
        assert_eq!(outer_map.get(&3), Some(&vec![1, 2, 5, ]));
        assert_eq!(outer_map.get(&4), Some(&vec![1, 2, 5, ]));
        assert_eq!(outer_map.get(&5), Some(&vec![1, 2, 3, 4]));
    }
    #[test]
    fn test_laplacian_smoothing() {
        let sqrt3  = 1.7320508075688772;
        let sqrt27 = sqrt3 * 3.0;

        let mut points = vec![
            Point { x: Vector3::new(6.0, 1.0, 0.0) }, // 0
            Point { x: Vector3::new(7.0, 0.0, 1.0) }, // 1
            Point { x: Vector3::new(8.0, -1.0, 0.0) }, // 2
            Point { x: Vector3::new(9.0, 0.0, -1.0) }, // 3
            Point { x: Vector3::new( 6.0, 0.0, 0.0) }, // 4
            Point { x: Vector3::new(-3.0,  sqrt27, 0.0) }, // 5
            Point { x: Vector3::new(-3.0, -sqrt27, 0.0) }, // 6
        ];
        let inner_index = vec![0, 1, 2, 3];
        let neighbor_map = {
            let mut neighbor_map: HashMap<usize, Vec<usize>> = HashMap::new();
            neighbor_map.insert(0, vec![1, 2, 3]);
            neighbor_map.insert(1, vec![0, 4, 5]);
            neighbor_map.insert(2, vec![0, 4, 6]);
            neighbor_map.insert(3, vec![0, 5, 6]);
            neighbor_map
        };
        let iteration: usize = 50;
        for _ in 0..iteration {
            points = laplacian_smoothing(points.clone(), inner_index.clone(), neighbor_map.clone());
        }
        assert_eq!(points[0], Point { x: Vector3::new(0.0, 0.0, 0.0) });
        assert_eq!(points[1], Point { x: Vector3::new(1.0, sqrt3, 0.0) });
        assert_eq!(points[2], Point { x: Vector3::new(1.0, -sqrt3, 0.0) });
        assert_eq!(points[3], Point { x: Vector3::new(-2.0, 0.0, 0.0) });
    }

    #[test]
    fn test_point_direction(){
        let point = Point { x: Vector3::new(10.0, 10.0, 10.0) };
        let direction: Vector3<f32> = point.direction(Vector3::new(6.0, 7.0, 10.0));
        assert_eq!(direction, Vector3::new(0.8, 0.6, 0.0));
    }
    #[test]
    fn test_point_orthogonal(){
        let point = Point { x: Vector3::new(10.0, 10.0, 10.0) };
        let center:Vector3<f32> = Vector3::new(6.0, 7.0, 10.0);
        let axis:Vector3<f32> = Vector3::new(1.0, 0.0, 0.0);
        let (new_center, direction) = point.orthogonal(center, axis);
        assert_eq!(new_center, Vector3::new(10.0, 7.0, 10.0));
        assert_eq!(direction, Vector3::new(0.0, 1.0, 0.0));
    }
    #[test]
    fn test_point_orthogonal_cone(){
        let sqrt1_3  = 1.7320508075688772 / 3.0;
        let point = Point { x: Vector3::new(10.0, 10.0, 10.0) };
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
        assert_eq!(new_point3, Point { x: Vector3::new(11.0, 10.0, 9.0) });
    }
    #[test]
    fn test_point_project_on_cylinder() {
        let point = Point { x: Vector3::new(-8.0, -6.0, 1.0) };
        let center:Vector3<f32> = Vector3::zeros();
        let axis:Vector3<f32> = Vector3::new(0.0, 0.0, 1.0);
        let radius = 5.0;
        let new_point = point.project_on_cylinder(center, axis, radius);
        assert_eq!(new_point, Point { x: Vector3::new(-4.0, -3.0, 1.0) });
    }
    #[test]
    fn test_point_project_on_sphire() {
        let point = Point { x: Vector3::new(-8.0, -6.0, 0.0) };
        let center:Vector3<f32> = Vector3::zeros();
        let radius = 5.0;
        let new_point = point.project_on_sphire(center, radius);
        assert_eq!(new_point, Point { x: Vector3::new(-4.0, -3.0, 0.0) });
    }
    #[test]
    fn test_laplacian_smoothing_with_axis_normalizing() {
        let sqrt3  = 1.7320508075688772;

        let mut points = vec![
            Point { x: Vector3::new(0.0, 2.0, 0.0) }, // 0
            Point { x: Vector3::new(0.0, 0.0, 1.0) }, // 1
            Point { x: Vector3::new(0.0, 0.0, 2.0) }, // 2
            Point { x: Vector3::new(2.0, 0.0, 0.0) }, // 3
        ];
        let inner_index = vec![1, 2, ];
        let outer_map = {
            let mut outer_map: HashMap<usize, Vec<usize>> = HashMap::new();
            outer_map.insert(1, vec![0, 2]);
            outer_map.insert(2, vec![1, 3]);
            outer_map
        };
        let iteration: usize = 50;
        for _ in 0..iteration {
            points = laplacian_smoothing_with_axis_normalizing(points.clone(), inner_index.clone(), outer_map.clone(), Vector3::zeros(), Vector3::new(0.0, 0.0, 1.0), 2.0);
        }

        assert_ne!(points[1], Point { x: Vector3::new(1.0, sqrt3, 0.0) });
        assert_ne!(points[2], Point { x: Vector3::new(sqrt3, 1.0, 0.0) });
    }

}



