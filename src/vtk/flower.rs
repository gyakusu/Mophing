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
}
