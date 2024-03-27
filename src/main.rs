use mophing::vtk::io;

fn main() {

    let old_file_path = "data/Tetra.vtu";
    let mut mesh = io::read_all(old_file_path);

    let mut twice_point = mesh.points.clone();
    twice_point.iter_mut().for_each(|p| {
        p.mul(2.0);
    });
    mesh.points.extend(twice_point);

    mesh.smooth_inner(100);

    let new_file_path = old_file_path.replace(".vtu", "_smoothed.vtu"); 
    io::copy_vtk_and_replace_point(old_file_path, &new_file_path, &mut mesh.points);

}




