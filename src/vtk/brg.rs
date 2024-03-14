use super::*;

pub struct Bearing {
}

impl Bearing {
    pub fn get_inner_points(&self, mesh: Mesh) -> Vec<Point> {
        let points = mesh.get_inner_points()
    }
}

