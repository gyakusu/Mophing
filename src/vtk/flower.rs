use super::face::Face;

// 面の2等分面を定義するインデックスを保持するクラス
#[derive(Debug, PartialEq, Clone, Copy, Hash, Eq)]
pub struct Flower {
    bottom: Face,       // 花の底面のインデックス
    petal: [usize; 3],  // 花びら3枚のインデックス
    is_front: bool,     // 花びらが表面か裏面か

}
impl Flower {
    pub fn new(bottom: Face, petal: [usize; 3], is_front: bool) -> Self {
        Self {
            bottom,
            petal,
            is_front,
        }
    }
    pub fn from_vec(index: [usize; 6]) -> Self {
        let (bottom, is_front) = Face::new([index[0], index[1], index[2]]);
        Self {
            bottom,
            petal: [index[3], index[4], index[5]],
            is_front,
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
        if self.is_front {
            [i0, i1, i2]
        } else {
            [i1, i0, i2]
        }
    }
    pub fn bottom(&self) -> [usize; 3] {
        self.bottom.as_vec()
    }
    pub fn flip(&self) -> f32 {
        if self.is_front { 1.0 } else { -1.0 }
    }
}
