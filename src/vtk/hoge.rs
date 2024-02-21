use super::reader;
use super::vtk;

pub fn fuga()-> Vec<vtk::Point>{
    let p = reader::read_point("data/Tetra.vtu");
    p
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hoge() {

        let points = fuga();
        println!("{:?}", points.first());
        println!("{:?}", points.last());
        println!("{:?}", points.len());

        assert!(true);
    }
}

