use super::point::Point;

pub struct Face {
    pub points: Vec<Point>,
    pub faces: Vec<i64>,
    pub faceoffsets: Vec<i64>,
}
