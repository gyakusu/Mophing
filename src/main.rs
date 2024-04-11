use mophing::vtk::io;
use mophing::vtk::brg::Brg;

fn main() {

    let origin_path = "data/Tetra.vtu";
    let index_path = "data/face_and_edge_index.xml";
    let write_path = "data/Tetra_linspace.vtu";
    // let setting_path = "data/Tetra_setting.xml";
    // let mesh = io::read_vtk_and_setting(origin_path, setting_path);

    // let origin_path = "data/Tetra_Cage.vtu";
    // let index_path = "data/face_and_edge_index_cage.xml";
    // let write_path = "data/Tetra_linspace_cage.vtu";
    // let setting_path = "data/Tetra_setting_cage.xml";

    let mesh = io::read_vtk(origin_path);

    let mut brg = Brg::sample(&mesh);

    let section0 = "edge";
    let section1 = "face";
    let edges = io::read_index_from_xml(index_path, section0).unwrap();
    let faces = io::read_index_from_xml(index_path, section1).unwrap();

    brg.set_edge_and_face(edges, faces);

    brg.linspace_all();

    brg.project_all();

    brg.smooth_face_all(20);
    brg.smooth_inner(10);
    brg.smooth_face_all(100);
    brg.smooth_inner(50);

    io::copy_vtk_and_replace_point(origin_path, write_path, &brg.get_points());

}




