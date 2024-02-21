#[derive(Debug, PartialEq)]
pub struct Point {
    pub x: [f32; 3],
}

#[derive(Debug, PartialEq)]
pub struct Tetra {
    pub index: [i64; 4],
}
