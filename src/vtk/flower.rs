use super::face::Face;

// 面の2等分面を定義するインデックスを保持するクラス
#[derive(Debug, PartialEq, Clone, Copy, Hash, Eq)]
pub struct Flower {
    bottom: Face, // 花の底面のインデックス
    petal: [usize; 3],  // 花びら3枚のインデックス
    
}
impl Flower {
    pub fn new(bottom: Face, petal: [usize; 3]) -> Self {
        Self {
            bottom,
            petal,
        }
    }
    pub fn from_vec(index: [usize; 6]) -> Self {
        Self {
            bottom: Face::new([index[0], index[1], index[2]]),
            petal: [index[3], index[4], index[5]],
        }
    }
    pub fn get(&self, i: usize) -> [usize; 3] {
        let j: [usize; 2] = match i {
            0 => [1, 2],
            1 => [2, 0],
            2 => [0, 1],
            _ => [3, 3],
        };
        let i0 = self.bottom.get(j[0]);
        let i1 = self.bottom.get(j[1]);
        let i2 = self.petal[i];
        [i0, i1, i2]
    }
    pub fn bottom(&self) -> [usize; 3] {
        self.bottom.as_vec()
    }
}
