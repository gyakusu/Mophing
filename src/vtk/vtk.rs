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
}

impl Mesh {
    pub fn new(points: Vec<Point>, tetras: Vec<Tetra>) -> Self {
        
        let inner_index = inner_point_indices(&tetras);

        Self {
            points,
            tetras,
            inner_index: inner_index,
        }
    }
}

pub fn inner_point_indices(tetras: &Vec<Tetra>) -> Vec<i64> {

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

    let all_indices: Vec<[i64; 3]> = faces.iter().map(|face| face.index).collect();
    
    let inner_indices: Vec<i64> = all_indices
    .into_iter()
    .flat_map(|index| index.iter())
    .filter(|&index| !outer_faces.iter().any(|face| face.index.contains(index)))
    .cloned()
    .collect();
    
    inner_indices
}

// pub fn found_neighbors();


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_inner_point_indices() {

        let tetras = vec![
            Tetra::new([0, 1, 2, 3]).unwrap(),
            Tetra::new([0, 1, 2, 4]).unwrap(),
            Tetra::new([0, 1, 3, 4]).unwrap(),
            Tetra::new([0, 2, 3, 4]).unwrap(),
        ];
        let inner_point_indices = inner_point_indices(&tetras);
        assert_eq!(inner_point_indices, vec![0, ]);
    }
}


