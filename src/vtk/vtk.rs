use std::collections::HashSet;
use std::collections::HashMap;

#[derive(Debug, PartialEq)]
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
            return Err("p_index must have a length of 4");
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
}

impl Mesh {
    pub fn new(points: Vec<Point>, tetras: Vec<Tetra>) -> Self {
        
        Self {
            points,
            tetras,
        }
    }

    pub fn outer_point_indices(&self) -> Vec<usize> {
        let mut faces: Vec<Face> = Vec::new();
        for tetra in &self.tetras {
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
        let mut face_counts: HashMap<Face, usize> = HashMap::new();
        for face in &faces {
            *face_counts.entry(face.clone()).or_insert(0) += 1;
        }
        let not_shared_faces: HashSet<Face> = face_counts.into_iter()
            .filter(|(_, count)| *count == 1)
            .map(|(face, _)| face)
            .collect();
        let mut indices: HashSet<usize> = HashSet::new();
        for face in not_shared_faces {
            indices.extend(face.index.iter().map(|&i| i as usize));
        }
        let mut indices_vec = Vec::from_iter(indices); 
        indices_vec.sort();
        indices_vec

    }

}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_outer_point_indices() {
        let points = vec![
            Point { x: [1.0, 1.0, 1.0] },
            Point { x: [2.0, 0.0, 0.0] },
            Point { x: [0.0, 2.0, 0.0] },
            Point { x: [0.0, 0.0, 2.0] },
            Point { x: [2.0, 2.0, 2.0] },
        ];
        let tetras = vec![
            Tetra::new([0, 1, 2, 3]).unwrap(),
            Tetra::new([0, 1, 2, 4]).unwrap(),
            Tetra::new([0, 1, 3, 4]).unwrap(),
            Tetra::new([0, 2, 3, 4]).unwrap(),
        ];
        let mesh = Mesh::new(points, tetras);
        let outer_point_indices = mesh.outer_point_indices();
        assert_eq!(outer_point_indices, vec![1, 2, 3, 4]);
    }
}








// #[derive(Debug, Eq, Hash, PartialEq, Copy, Clone, PartialOrd, Ord)]
// pub struct Node {
//     pub p_index: [i64; 2],
// }

// let mut nodes = Vec::new();

// for tetra in &tetras {
//     for i in 0..4 {
//         for j in (i + 1)..4 {
//             let (min_index, max_index) = if tetra.p_index[i] < tetra.p_index[j] {
//                 (tetra.p_index[i], tetra.p_index[j])
//             } else {
//                 (tetra.p_index[j], tetra.p_index[i])
//             };
//             let p_index = [min_index, max_index];
//             let node = Node { p_index };
//             nodes.push(node);
//         }
//     }
// }

// let mut node_map: HashMap<Node, Vec<usize>> = HashMap::new();

// for (i, node) in nodes.iter().enumerate() {
//     node_map.entry(*node).or_default().push(i);
// }

// let mut removed_indices = Vec::new();
// let mut unique_nodes = Vec::new();

// for (node, mut indices) in node_map {
//     indices.sort();
//     let min_index = indices[0];
//     removed_indices.extend(indices.into_iter().skip(1));
//     unique_nodes.push((min_index, node));
// }

// // unique_nodes.sort_by_key(|&(i, _)| i);
// unique_nodes.sort_by_key(|&(_,i)| i);

// let original_indices: Vec<usize> = unique_nodes.iter().map(|&(i, _)| i).collect();

// let mut indexed_nodes: Vec<(usize, &Node)> = nodes.iter().enumerate().collect();

// indexed_nodes.sort_by(|a, b| {
//     let cmp = a.1.p_index[0].cmp(&b.1.p_index[0]);
//     if cmp == std::cmp::Ordering::Equal {
//         a.1.p_index[1].cmp(&b.1.p_index[1])
//     } else {
//         cmp
//     }
// });

// indexed_nodes.dedup_by(|a, b| a.1.p_index == b.1.p_index);

// let original_indices: Vec<usize> = indexed_nodes.iter().map(|(i, _)| *i).collect();
