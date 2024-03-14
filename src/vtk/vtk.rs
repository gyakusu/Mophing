use std::collections::HashSet;
use std::collections::HashMap;

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Point {
    pub x: [f32; 3],
}

#[derive(Debug, Eq, Hash, PartialEq, Clone)]
pub struct Face {
    pub index: [i64; 3],
}

#[derive(Debug, PartialEq)]
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
    pub inner_index: Vec<i64>,
    pub neighbor_index: HashMap<i64, Vec<i64>>,
}

impl Mesh {
    pub fn new(points: Vec<Point>, tetras: Vec<Tetra>) -> Self {
        
        let inner_index = find_inners(&tetras);
        let neighbor_index = find_neighbors(&tetras);

        Self {
            points,
            tetras,
            inner_index: inner_index,
            neighbor_index: neighbor_index,
        }
    }
    pub fn laplacian_smoothing(&self, iteration: i64) {
        // let new_points = laplacian_smoothing(self.points, self.inner_index, self.neighbor_index, iteration);
        // self.points = new_points;
    }
}

pub fn find_inners(tetras: &Vec<Tetra>) -> Vec<i64> {

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
    let mut face_counts: HashMap<Face, i64> = HashMap::new();
    for face in &faces {
        *face_counts.entry(face.clone()).or_insert(0) += 1;
    }
    let outer_faces: HashSet<Face> = face_counts.into_iter()
        .filter(|(_, count)| *count == 1)
        .map(|(face, _)| face)
        .collect();

    let outer_indices: HashSet<i64> = outer_faces.iter()
        .flat_map(|face| face.index.iter())
        .cloned()
        .collect();

    let all_indices: HashSet<i64> = faces.iter()
        .flat_map(|face| face.index.iter())
        .cloned()
        .collect();

    let find_inners: Vec<i64> = all_indices.difference(&outer_indices).cloned().collect();
    find_inners
}

pub fn find_neighbors(tetras: &Vec<Tetra>) -> HashMap<i64, Vec<i64>> {
    let mut neighbor_index: HashMap<i64, Vec<i64>> = HashMap::new();
    for tetra in tetras {
        for i in 0..4 {
            for j in 0..4 {
                if i != j {
                    neighbor_index.entry(tetra.index[i]).or_insert_with(Vec::new).push(tetra.index[j]);
                }
            }
        }
    }
    for (_key, value) in neighbor_index.iter_mut() {
        value.sort_unstable(); // ソート
        value.dedup(); // 重複の削除
    }
    neighbor_index
}

pub fn laplacian_smoothing(points: Vec<Point>, inner_index: Vec<i64>, neighbor_index: HashMap<i64, Vec<i64>>, iteration: i64) -> Vec<Point> {
    let mut new_points = points.clone();

    for _ in 0..iteration {
        for &i in &inner_index {
            let mut sum = [0.0, 0.0, 0.0];
            let neighbors = &neighbor_index[&i];
            for j in neighbors {
                sum[0] += new_points[*j as usize].x[0];
                sum[1] += new_points[*j as usize].x[1];
                sum[2] += new_points[*j as usize].x[2];
            }
            let n = neighbors.len() as f32;
            new_points[i as usize].x[0] = sum[0] / n;
            new_points[i as usize].x[1] = sum[1] / n;
            new_points[i as usize].x[2] = sum[2] / n;
        }
    }
    new_points
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_find_inners() {

        let tetras = vec![
            Tetra::new([0, 1, 2, 3]).unwrap(),
            Tetra::new([0, 1, 2, 4]).unwrap(),
            Tetra::new([0, 1, 3, 4]).unwrap(),
            Tetra::new([0, 2, 3, 4]).unwrap(),
        ];
        let find_inners = find_inners(&tetras);
        assert_eq!(find_inners, vec![0, ]);
    }
    #[test]
    fn test_find_neighbors() {
        let tetras = vec![
            Tetra::new([0, 1, 3, 5]).unwrap(),
            Tetra::new([0, 1, 3, 6]).unwrap(),
            Tetra::new([0, 1, 4, 5]).unwrap(),
            Tetra::new([0, 1, 4, 6]).unwrap(),
            Tetra::new([0, 2, 3, 5]).unwrap(),
            Tetra::new([0, 2, 3, 6]).unwrap(),
            Tetra::new([0, 2, 4, 5]).unwrap(),
            Tetra::new([0, 2, 4, 6]).unwrap(),
        ];
        let neighbor_index = find_neighbors(&tetras);
        assert_eq!(neighbor_index.get(&0), Some(&vec![1, 2, 3, 4, 5, 6]));
        assert_eq!(neighbor_index.get(&1), Some(&vec![0, 3, 4, 5, 6]));
        assert_eq!(neighbor_index.get(&2), Some(&vec![0, 3, 4, 5, 6]));
        assert_eq!(neighbor_index.get(&3), Some(&vec![0, 1, 2, 5, 6]));
        assert_eq!(neighbor_index.get(&4), Some(&vec![0, 1, 2, 5, 6]));
        assert_eq!(neighbor_index.get(&5), Some(&vec![0, 1, 2, 3, 4]));
        assert_eq!(neighbor_index.get(&6), Some(&vec![0, 1, 2, 3, 4]));
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
        let neighbor_index = {
            let mut neighbor_index: HashMap<i64, Vec<i64>> = HashMap::new();
            neighbor_index.insert(0, vec![1, 2, 3]);
            neighbor_index.insert(1, vec![0, 4, 5]);
            neighbor_index.insert(2, vec![0, 4, 6]);
            neighbor_index.insert(3, vec![0, 5, 6]);
            neighbor_index
        };
        let new_points = laplacian_smoothing(points, inner_index, neighbor_index, 100);
        assert_eq!(new_points[0], Point { x: [ 0.0,    0.0, 0.0] });
        assert_eq!(new_points[1], Point { x: [ 1.0,  sqrt3, 0.0] });
        assert_eq!(new_points[2], Point { x: [ 1.0, -sqrt3, 0.0] });
        assert_eq!(new_points[3], Point { x: [-2.0,    0.0, 0.0] });
    }
}


