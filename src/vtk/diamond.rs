// 面の2等分面を定義するインデックスを保持するクラス
#[derive(Debug, PartialEq, Clone, Copy, Hash, Eq)]
pub struct Diamond {
    pub minor: [usize; 2], // 全てのインデックスと繋がるインデックス
    pub major: [usize; 2], // 残りの2つのインデックス
}
impl Diamond {
    pub fn new(minor0: usize, minor1: usize, major0: usize, major1: usize) -> Self {
        Self {
            minor: [minor0, minor1],
            major: [major0, major1],
        }
    }
}
