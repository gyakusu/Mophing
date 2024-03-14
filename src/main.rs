use mophing::vtk::reader;
use mophing::vtk::vtk;
use mophing::vtk::vtk::Mesh;


fn main() {

    let file_path = "data/Tetra.vtu";
    let mesh = reader::read_all(file_path);

    mesh.laplacian_smoothing(10);

    let new_file_path = file_path.replace(".vtu", "_smoothed.vtu"); 

    println!("{:?}", mesh.points.len());
}




