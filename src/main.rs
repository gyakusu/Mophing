use mophing::vtk::reader;
// use mophing::vtk::vtk::Mesh;


fn main() {
    let file_path = "data/Tetra.vtu";
    let mesh = reader::read_all(file_path);
    let outer_point_indices = mesh.outer_point_indices();
    println!("{:?}", mesh.points.len());
    println!("{:?}", outer_point_indices.len());
}




