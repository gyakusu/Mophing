use mophing::vtk::io;

fn main() {

    let old_file_path = "data/Tetra.vtu";
    let mut mesh = io::read_all(old_file_path);

    let mut twice_point = mesh.get_inner_points();
    twice_point.iter_mut().for_each(|p| {p.x[0] *= 2.0; p.x[1] *= 2.0; p.x[2] *= 2.0;});
    mesh.set_inner_points(twice_point);

    mesh.laplacian_smoothing(100);

    let new_file_path = old_file_path.replace(".vtu", "_smoothed.vtu"); 
    io::copy_vtk_and_replace_point(old_file_path, &new_file_path, &mut mesh.points);

}




