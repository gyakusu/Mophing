use super::*;
use super::point::Point;

pub struct LineEdge {
    index: Vec<i64>,
}
pub struct CurveEdge {
    index: Vec<i64>,
    axis: [f32; 3],
    radius: f32,
}
pub struct PocketEdge {
    index: Vec<i64>,
    center: [f32; 3],
    radius: f32,
}
pub trait Edge {
    fn new(index: Vec<i64>) -> Self;
    fn as_vec(&self) -> Vec<i64>;
    fn linspace(&self, points: &Vec<Point>, start: f32, end: f32) -> Vec<Point>;
}
impl Edge for PocketEdge {
    fn new(index: Vec<i64>) -> Self {
        Self {
            index,
        }
    }
    fn as_vec(&self) -> Vec<i64> {
        self.index.clone()
    }
}

pub struct CurveFace

pub struct Cage {
    pub mesh: Mesh,
    pub in_pocket: PocketEdge,
    pub out_pocket: PocketEdge,
    pub in_line: Edge,
    pub out_line: Edge,
    pub face_pocket: Edge,
    

}

impl Cage {

}

