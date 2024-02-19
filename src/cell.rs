use crate::point::Point;

pub struct Cell {
    pub points: Vec<Point>,
    pub connectivity: Vec<i64>,
    pub offsets: Vec<i64>,
}

// face.rs

