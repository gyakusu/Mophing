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
    pub fn as_i64(&self) -> [usize; 3] {
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
        assert_eq!(face.as_i64(), [1, 2, 3]);
    }

}


