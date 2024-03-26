use std::collections::HashSet;
use std::collections::HashMap;

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Point {
    x: [f32; 3],
}

impl Point {
    pub fn new(x: [f32; 3]) -> Self {
        Self {
            x: x,
        }
    }
    pub fn on_circle(center: [f32;3], direction: [f32;3], radius: f32) -> Self {
        Self {
            x: [center[0] + direction[0] * radius, center[1] + direction[1] * radius, center[2] + direction[2] * radius]
        }
    }
    pub fn as_f32(&self) -> [f32;3] {
        self.x
    }
    pub fn add(&self, a: [f32;3]) -> Self {
        Self {
            x: [self.x[0] + a[0], self.x[1] + a[1], self.x[2] + a[2]]
        }
    }
    pub fn mul(&self, a: f32) -> Self {
        Self {
            x: [self.x[0] * a, self.x[1] * a, self.x[2] * a],
        }
    }
    pub fn direction(&self, center: [f32;3]) -> [f32;3] {
        let dx = self.x[0] - center[0];
        let dy = self.x[1] - center[1];
        let dz = self.x[2] - center[2];
        let distance = (dx * dx + dy * dy + dz * dz).sqrt();
        if distance != 0.0 {[dx / distance, dy / distance, dz / distance]} else { [0.0, 0.0, 0.0] }
    }
    pub fn orthogonal(&self, center: [f32;3], axis: [f32;3]) -> ([f32;3], [f32;3]) {
        let dx = self.x[0] - center[0];
        let dy = self.x[1] - center[1];
        let dz = self.x[2] - center[2];
        let dot = dx * axis[0] + dy * axis[1] + dz * axis[2];
        let new_center = [center[0] + axis[0] * dot, center[1] + axis[1] * dot, center[2] + axis[2] * dot];
        let orth = [dx - dot * axis[0], dy - dot * axis[1], dz - dot * axis[2]];
        let distance = (orth[0] * orth[0] + orth[1] * orth[1] + orth[2] * orth[2]).sqrt();
        (new_center, if distance != 0.0 {[orth[0] / distance, orth[1] / distance, orth[2] / distance]} else { [0.0, 0.0, 0.0] })
    }
}

#[derive(Debug, Eq, Hash, PartialEq, Clone, Copy)]
pub struct Face {
    pub index: [i64; 3],
}

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Tetra {
    pub index: [i64; 4],
}

impl Tetra {
    pub fn new(index: [i64; 4]) -> Result<Self, &'static str> {
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

#[derive(Debug, PartialEq)]
pub struct Mesh {
    pub points: Vec<Point>,
    pub tetras: Vec<Tetra>,
    pub faces: Vec<Face>,
    pub outer_index: Vec<i64>,
    pub inner_index: Vec<i64>,
    pub neighbor_map: HashMap<i64, Vec<i64>>,
    pub outer_map: HashMap<i64, Vec<i64>>,
}

impl Mesh {
    pub fn new(points: Vec<Point>, tetras: Vec<Tetra>) -> Self {
        
        let faces: Vec<Face> = tetras_to_faces(&tetras);   
        let outer_index: Vec<i64> = find_outer_index(faces.clone());
        let inner_index: Vec<i64> = find_inner_index(faces.clone(), outer_index.clone());
        let neighbor_map: HashMap<i64, Vec<i64>> = find_neighbors(&tetras);
        let outer_map: HashMap<i64, Vec<i64>> = find_outer_neighbors(&neighbor_map, &inner_index);

        Self {
            points: points,
            tetras: tetras,
            faces: faces,
            outer_index: outer_index,
            inner_index: inner_index,
            neighbor_map: neighbor_map,
            outer_map: outer_map,
        }
    }
    pub fn smooth_inner(&mut self, iteration: i64) {
        let new_points = laplacian_smoothing(self.points.clone(), self.inner_index.clone(), self.neighbor_map.clone(), iteration);
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
pub fn find_outer_index(faces: Vec<Face>) -> Vec<i64> {

    let mut face_counts: HashMap<Face, i64> = HashMap::new();
    for face in &faces {
        *face_counts.entry(face.clone()).or_insert(0) += 1;
    }
    let shared_faces: HashSet<Face> = face_counts.iter()
        .filter(|&(_face, &count)| count == 2)
        .map(|(face, _count)| face.clone())
        .collect();
    let not_shared_faces: Vec<Face> = faces.iter()
        .filter(|face| !shared_faces.contains(face))
        .cloned()
        .collect();
    let outer_map: HashSet<i64> = not_shared_faces.iter()
        .flat_map(|face| face.index.iter())
        .cloned()
        .collect();
    let mut outer_index: Vec<i64> = outer_map.iter().cloned().collect();
    outer_index.sort();
    outer_index
}

pub fn find_inner_index(faces: Vec<Face>, outer_index: Vec<i64>) -> Vec<i64> {

    let all_index: HashSet<i64> = faces.iter()
        .flat_map(|face| face.index.iter())
        .cloned()
        .collect();
    let mut inner_index: Vec<i64> = all_index.difference(&outer_index.iter().cloned().collect())
        .cloned()
        .collect();
    inner_index.sort();
    inner_index
}

pub fn find_neighbors(tetras: &Vec<Tetra>) -> HashMap<i64, Vec<i64>> {
    let mut neighbor_map: HashMap<i64, Vec<i64>> = HashMap::new();
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

pub fn find_outer_neighbors(neighbor_map: &HashMap<i64, Vec<i64>>, inner_index: &[i64]) -> HashMap<i64, Vec<i64>> {
    let mut outer_map: HashMap<i64, Vec<i64>> = HashMap::new();
    for (key, value) in neighbor_map {
        let filtered_value: Vec<i64> = value.iter().filter(|&index| !inner_index.contains(index)).cloned().collect();
        outer_map.insert(*key, filtered_value);
    }
    outer_map
}

pub fn laplacian_smoothing(points: Vec<Point>, inner_index: Vec<i64>, neighbor_map: HashMap<i64, Vec<i64>>, iteration: i64) -> Vec<Point> {
    let mut new_points = points.clone();

    for _ in 0..iteration {
        for &i in &inner_index {
            let neighbors = &neighbor_map[&i];
            let sum = neighbors.iter().fold(Point { x: [0.0, 0.0, 0.0] }, |sum, &j| sum.add(new_points[j as usize].as_f32()));
            new_points[i as usize] = sum.mul(1.0 / neighbors.len() as f32);
        }
    }
    new_points
}

pub fn laplacian_smoothing_with_center_normalizing(points: Vec<Point>, inner_index: Vec<i64>, neighbor_map: HashMap<i64, Vec<i64>>, center: [f32; 3], radius: f32, iteration: i64) -> Vec<Point> {
    let mut new_points = points.clone();

    for _ in 0..iteration {
        for &i in &inner_index {
            let neighbors = &neighbor_map[&i];
            let sum = neighbors.iter().fold(Point { x: [0.0, 0.0, 0.0] }, |sum, &j| sum.add(new_points[j as usize].as_f32()));
            let direction = sum.direction(center);
            new_points[i as usize] = Point::on_circle(center, direction, radius);
        }
    }
    new_points
}
pub fn laplacian_smoothing_with_axis_normalizing(points: Vec<Point>, inner_index: Vec<i64>, neighbor_map: HashMap<i64, Vec<i64>>, center: [f32; 3], axis: [f32; 3], radius: f32, iteration: i64) -> Vec<Point> {
    let mut new_points = points.clone();

    for _ in 0..iteration {
        for &i in &inner_index {
            let neighbors = &neighbor_map[&i];
            let sum = neighbors.iter().fold(Point { x: [0.0, 0.0, 0.0] }, |sum, &j| sum.add(new_points[j as usize].as_f32()));
            let mean = sum.mul(1.0 / neighbors.len() as f32);
            let (new_center, direction) = mean.orthogonal(center, axis);
            new_points[i as usize] = Point::on_circle(new_center, direction, radius);
        }
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
        let faces: Vec<Face> = tetras_to_faces(&tetras);
        let outer_index = find_outer_index(faces);
        assert_eq!(outer_index, vec![1, 2, 3, 4]);
    }

    #[test]
    fn test_find_inner_index() {
        let tetras = TETRAS0.to_vec();
        let faces: Vec<Face> = tetras_to_faces(&tetras);
        let outer_index: Vec<i64> = vec![1, 2, 3, 4];
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

        let points = vec![
            Point { x: [6.0, 1.0, 0.0] }, // 0
            Point { x: [7.0, 0.0, 1.0] }, // 1
            Point { x: [8.0, -1.0, 0.0] }, // 2
            Point { x: [9.0, 0.0, -1.0] }, // 3
            Point { x: [ 6.0, 0.0, 0.0] }, // 4
            Point { x: [-3.0,  sqrt27, 0.0] }, // 5
            Point { x: [-3.0, -sqrt27, 0.0] }, // 6
        ];
        let inner_index = vec![0, 1, 2, 3];
        let neighbor_map = {
            let mut neighbor_map: HashMap<i64, Vec<i64>> = HashMap::new();
            neighbor_map.insert(0, vec![1, 2, 3]);
            neighbor_map.insert(1, vec![0, 4, 5]);
            neighbor_map.insert(2, vec![0, 4, 6]);
            neighbor_map.insert(3, vec![0, 5, 6]);
            neighbor_map
        };
        let new_points = laplacian_smoothing(points, inner_index, neighbor_map, 100);

        assert_eq!(new_points[0], Point { x: [ 0.0,    0.0, 0.0] });
        assert_eq!(new_points[1], Point { x: [ 1.0,  sqrt3, 0.0] });
        assert_eq!(new_points[2], Point { x: [ 1.0, -sqrt3, 0.0] });
        assert_eq!(new_points[3], Point { x: [-2.0,    0.0, 0.0] });
    }

    #[test]
    fn test_point_direction(){
        let point = Point { x: [10.0, 10.0, 10.0] };
        let direction = point.direction([6.0, 7.0, 10.0]);
        assert_eq!(direction, [0.8, 0.6, 0.0]);
    }
    #[test]
    fn test_point_orthogonal(){
        let point = Point { x: [10.0, 10.0, 10.0] };
        let (new_center, direction) = point.orthogonal([6.0, 7.0, 10.0], [1.0, 0.0, 0.0]);
        assert_eq!(new_center, [10.0, 7.0, 10.0]);
        assert_eq!(direction, [0.0, 1.0, 0.0]);
    }
}


