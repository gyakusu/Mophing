use super::math::sort_three_numbers as sort_three_numbers;

#[derive(Debug, Eq, Hash, PartialEq, Clone, Copy)]
pub struct Face {
    index: [usize; 3],
}
impl Face {
    pub fn new(index: [usize; 3]) -> Self {
        #[cfg(debug_assertions)] {
            if index.len() != 3 {
                panic!("index must have a length of 3");
            }
        }
        let sorted_index = sort_three_numbers(index);
        Self {
            index: sorted_index,
        }
    }
    pub fn as_vec(&self) -> [usize; 3] {
        self.index
    }
    pub fn get(&self, i: usize) -> usize {
        #[cfg(debug_assertions)] {
            if i > 2 {
                panic!("Index out of bounds: {}", i);
            }
        }
        self.index[i]
    }
    pub fn neighbor(&self, new: usize, i: usize) -> (usize, Self) {
        #[cfg(debug_assertions)] {
            if i > 2 {
                panic!("Index out of bounds: {}", i);
            }
        }
        let remain = self.index[i];
        let j: [usize; 2] = match i {
            0 => [1, 2],
            1 => [0, 2],
            2 => [0, 1],
            _ => [3, 3],
        };
        let face = Face::new([self.index[j[0]], self.index[j[1]], new]);
        (remain, face)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_face_new() {
        let face = Face::new([1, 2, 3]);
        assert_eq!(face.get(0), 1);
        assert_eq!(face.get(1), 2);
        assert_eq!(face.get(2), 3);
    }
    #[test]
    fn test_face_as_i64() {
        let face = Face::new([1, 2, 3]);
        assert_eq!(face.as_vec(), [1, 2, 3]);
    }
    #[test]
    fn test_face_neighbor() {
        let face = Face::new([10, 20, 40]);
        let (remain0, new_face0) = face.neighbor(30, 0);
        assert_eq!(remain0, 10);
        assert_eq!(new_face0.as_vec(), [20, 30, 40]);
        let (remain1, new_face1) = face.neighbor(30, 1);
        assert_eq!(remain1, 20);
        assert_eq!(new_face1.as_vec(), [10, 30, 40]);
        let (remain2, new_face2) = face.neighbor(30, 2);
        assert_eq!(remain2, 40);
        assert_eq!(new_face2.as_vec(), [10, 20, 30]);
    }
}


