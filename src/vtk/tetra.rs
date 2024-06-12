use core::panic;
use std::collections::HashMap;
use std::collections::HashSet;
extern crate nalgebra as na;

use super::math::has_duplication_vec4;
use super::face::Face;
use super::flower::Flower;

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Tetra {
    pub index: [usize; 4],
}
impl Tetra {
    pub fn new(index: [usize; 4]) -> Result<Self, &'static str> {
        if index.len() != 4 {
            return Err("index must have a length of 4");
        }
        let has_duplication = has_duplication_vec4(&index);
        if has_duplication {
            return Err("index must not have duplication");
        }
        Ok(Self {
            index,
        })
    }
    pub fn remain_face(&self, i: usize) -> (usize, Face, bool) {
        #[cfg(debug_assertions)] {
            if i > 3 {
                panic!("Index out of bounds: {}", i);
            }
        }
        let j: [usize; 3] = match i {
            0 => [1, 2, 3],
            1 => [2, 0, 3],
            2 => [3, 0, 1],
            3 => [0, 2, 1],
            _ => [4, 4, 4],
        };
        let (face, is_front) = Face::new([self.index[j[0]], self.index[j[1]], self.index[j[2]]]);
        (self.index[i], face, is_front)
    }
    pub fn contain_and_remain(&self, face: &Face) -> (bool, usize) {
        let mut contain: [bool; 4] = [false; 4];

        for i in 0..4 {
            contain[i] = face.as_vec().contains(&self.index[i]);
        }
        let sum = contain.iter().fold(0, |sum, &b| sum + b as usize);
        let remain = contain.iter().position(|&b| !b).unwrap();
        (sum == 3, self.index[remain])
    }
}

pub fn tetras_to_face_map(tetras: &Vec<Tetra>) -> HashMap<Face, HashSet<(usize, bool)>> {
    let mut face_map: HashMap<Face, HashSet<(usize, bool)>> = HashMap::new();

    for tetra in tetras {
        for i in 0..4 {
            let (remain, face, is_front) = tetra.remain_face(i);
            if face_map.contains_key(&face) {
                face_map.get_mut(&face).unwrap().insert((remain, is_front));
            } else {
                face_map.insert(face, [(remain, is_front)].iter().cloned().collect());
            }
        }
    }
    face_map
}
pub fn find_surface_faces(face_map: &HashMap<Face, HashSet<(usize, bool)>>) -> HashSet<Face> {
    let mut surface_faces: HashSet<Face> = HashSet::new();
    for (face, remains) in face_map {
        if remains.len() == 1 {
            surface_faces.insert(*face);
        }
    }
    surface_faces
}
pub fn find_inner_index(face_map: &HashMap<Face, HashSet<(usize, bool)>>, surface_faces: &HashSet<Face>) -> HashSet<usize> {

    let mut inner_index: HashSet<usize> = HashSet::new();
    for (face, _) in face_map {
        inner_index.insert(face.get(0));
        inner_index.insert(face.get(1));
        inner_index.insert(face.get(2));
    }
    for surface_face in surface_faces {
        inner_index.remove(&surface_face.get(0));
        inner_index.remove(&surface_face.get(1));
        inner_index.remove(&surface_face.get(2));
    }
    inner_index.into_iter().collect()
 }

pub fn get_neighbor_map(tetras: &Vec<Tetra>) -> HashMap<usize, HashSet<usize>> {
    let mut neighbor_map: HashMap<usize, HashSet<usize>> = HashMap::new();
    for tetra in tetras {
        for i in 0..4 {
            for j in 0..4 {
                if i != j {
                    neighbor_map.entry(tetra.index[i]).or_insert_with(HashSet::new).insert(tetra.index[j]);
                }
            }
        }
    }
    neighbor_map
}
pub fn get_surface_map(surface_faces: &HashSet<Face>) -> HashMap<usize, HashSet<usize>> {
    let mut surface_map: HashMap<usize, HashSet<usize>> = HashMap::new();
    for face in surface_faces {
        for &i in &face.as_vec() {
            let neighbors: HashSet<usize> = face.as_vec().iter().cloned().filter(|&j| j != i).collect();
            surface_map.entry(i).or_insert_with(HashSet::new).extend(neighbors);
        }
    }
    surface_map
}
pub fn get_neighbors(neighbor_map: &HashMap<usize, HashSet<usize>>, i: usize) -> &HashSet<usize> {
    #[cfg(debug_assertions)] {
        match neighbor_map.get(&i) {
            Some(neighbors) => neighbors,
            None => {
                println!("Key not found in neighbor_map: {}", i);
                panic!("Key not found");
            }
        }
    }
    #[cfg(not(debug_assertions))] {
        neighbor_map.get(&i).unwrap()
    }
}
fn get_another(remains: &HashSet<(usize, bool)>, r0: usize) -> usize {
    if remains.len() != 2 {
        panic!("Remains must have a length of 2");
    }
    let r0_count = remains.iter().filter(|&x| x.0 == r0).count();
    if r0_count == 0 {
        panic!("Remains must contain r0");
    }
    if r0_count == 2 {
        panic!("Remains must contain only one r0");
    }
    let num_front = remains.iter().filter(|&x| x.1).count();
    if num_front != 1 {
        panic!("Remains must contain same front. num_front: {}", num_front);
    }
    for r in remains {
        if r.0 != r0 {
            return r.0;
        }
    }
    panic!("Key not found");
}
pub fn make_inverse_map(tetras: &Vec<Tetra>, inner_index: &HashSet<usize>, face_map: &HashMap<Face, HashSet<(usize, bool)>>) -> HashMap<usize, HashSet<Flower>> {
    let mut inverse_map: HashMap<usize, HashSet<Flower>> = HashMap::new();

    for tetra in tetras {
        for i in 0..4 {
            let (remain, face, is_front) = tetra.remain_face(i);
            if !inner_index.contains(&remain) {
                continue;
            }
            let mut petal: [usize; 3] = [0; 3];
            for j in 0..3 {
                let (r0, f0) = face.neighbor(remain, j);
                let remains = face_map.get(&f0).unwrap();
                let another = get_another(remains, r0);
                petal[j] = another;
            }
            let flower = Flower::new(face, petal, is_front);
            inverse_map.entry(remain).or_insert_with(HashSet::new).insert(flower);
        }
    }
    inverse_map
}

#[cfg(test)]
mod tests {
    use super::*;

    const TETRAS0 : [Tetra; 4] = [
        Tetra { index: [0, 1, 2, 3] },
        Tetra { index: [0, 1, 3, 4] },
        Tetra { index: [0, 1, 4, 2] },
        Tetra { index: [0, 2, 4, 3] },
    ];
    const TETRAS1 : [Tetra; 8] = [
        Tetra { index: [0, 1, 3, 5] },
        Tetra { index: [0, 3, 2, 5] },
        Tetra { index: [0, 2, 4, 5] },
        Tetra { index: [0, 4, 1, 5] },
        Tetra { index: [0, 1, 4, 6] },
        Tetra { index: [0, 4, 2, 6] },
        Tetra { index: [0, 2, 3, 6] },
        Tetra { index: [0, 3, 1, 6] },
    ];
    #[test]
    fn test_contain_and_remain() {
        let tetra = Tetra::new([0, 10, 20, 30]).unwrap();
        let (face0, _) = Face::new([0, 10, 20]);
        let (contain0, remain0) = tetra.contain_and_remain(&face0);
        assert_eq!(contain0, true);
        assert_eq!(remain0, 30);

        let (face1, _) = Face::new([0, 10, 40]);
        let (contain1, remain1) = tetra.contain_and_remain(&face1);
        assert_eq!(contain1, false);
        assert!(remain1 > usize::MIN);
    }
    #[test]
    fn test_tetras_to_face_map() {
        let tetras = TETRAS0.to_vec();
        let face_map = tetras_to_face_map(&tetras);
        let remain0 = face_map.get(&Face::new([0, 1, 2]).0).unwrap();
        assert!(remain0.contains(&(3, false)));
        assert!(remain0.contains(&(4, true)));

        let remain1 = face_map.get(&Face::new([1, 2, 3]).0).unwrap();
        assert!(remain1.contains(&(0, true)));
    }
    #[test]
    fn test_surface_faces() {
        let tetras = TETRAS0.to_vec();
        let face_map = tetras_to_face_map(&tetras);
        let surface_faces = find_surface_faces(&face_map);
        assert!(surface_faces.contains(&Face::new([1, 2, 3]).0));
        assert!(surface_faces.contains(&Face::new([1, 2, 4]).0));
        assert!(surface_faces.contains(&Face::new([1, 3, 4]).0));
        assert!(surface_faces.contains(&Face::new([2, 3, 4]).0));
    }
    #[test]
    fn test_find_inner_index() {
        let tetras = TETRAS0.to_vec();
        let face_map = tetras_to_face_map(&tetras);
        let surface_faces = find_surface_faces(&face_map);
        let inner_index = find_inner_index(&face_map, &surface_faces);
        assert!(inner_index.contains(&0));
    }
    #[test]
    fn test_find_neighbors() {
        let tetras = TETRAS1.to_vec();
        let neighbor_map = get_neighbor_map(&tetras);
        for &i in &[1, 2, 3, 4, 5, 6] {
            assert!(neighbor_map.get(&0).unwrap().contains(&i));
        }
        for &i in &[0, 3, 4, 5, 6] {
            assert!(neighbor_map.get(&1).unwrap().contains(&i));
            assert!(neighbor_map.get(&2).unwrap().contains(&i));
        }
        for &i in &[0, 1, 2, 5, 6] {
            assert!(neighbor_map.get(&3).unwrap().contains(&i));
            assert!(neighbor_map.get(&4).unwrap().contains(&i));
        }
        for &i in &[0, 1, 2, 3, 4] {
            assert!(neighbor_map.get(&5).unwrap().contains(&i));
            assert!(neighbor_map.get(&6).unwrap().contains(&i));
        }
    }
    #[test]
    fn test_find_surface_neighbors() {
        let tetras = TETRAS1.to_vec();
        let face_map = tetras_to_face_map(&tetras);
        let surface_faces = find_surface_faces(&face_map);
        let surface_map = get_surface_map(&surface_faces);
        for &i in &[3, 4, 5, 6] {
            assert!(surface_map.get(&1).unwrap().contains(&i));
            assert!(surface_map.get(&2).unwrap().contains(&i));
        }
        for &i in &[1, 2, 5, 6] {
            assert!(surface_map.get(&3).unwrap().contains(&i));
            assert!(surface_map.get(&4).unwrap().contains(&i));
        }
        for &i in &[1, 2, 3, 4] {
            assert!(surface_map.get(&5).unwrap().contains(&i));
            assert!(surface_map.get(&6).unwrap().contains(&i));
        }
    }
    #[test]
    fn test_make_inverse_map() {
        let tetras = TETRAS1.to_vec();
        let face_map = tetras_to_face_map(&tetras);
        let surface_faces = find_surface_faces(&face_map);
        let inner_index = find_inner_index(&face_map, &surface_faces);
        let inverse_map = make_inverse_map(&tetras, &inner_index, &face_map);

        let d0_len = inverse_map.get(&0).unwrap_or(&HashSet::new()).len();
        assert_eq!(d0_len, 8);

        for h in inverse_map.get(&0).unwrap() {
            println!("{:?}", h);
        }

        let d1_len = inverse_map.get(&1).unwrap_or(&HashSet::new()).len();
        assert_eq!(d1_len, 0);
    }
}

