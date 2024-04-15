use nalgebra as na;
use na::Vector3;
use std::f32::consts::PI;

use mophing::vtk::io;
use mophing::vtk::brg::Brg;
use mophing::vtk::brg::{CageParameter, PocketParameter, NeckParameter};

fn main() {

    let origin_path = "data/Tetra.vtu";
    let index_path = "data/face_and_edge_index.xml";
    let write_path = "data/Tetra_linspace.vtu";
    // let setting_path = "data/Tetra_setting.xml";
    // let mesh = io::read_vtk_and_setting(origin_path, setting_path);

    let mesh = io::read_vtk(origin_path);

    let param = CageParameter {
        axis: Vector3::new(0.0, 0.0, 1.0),
        theta0: -PI / 6.0,
        theta1:  PI / 6.0,
        r0: 2.345e-3 * 1.05,
        r1: 2.850e-3 * 0.95,
        h0: 0.93e-3 * 0.95,
        h1: 2.10e-3 * 1.05,
        bevel: 0.10e-3 * 0.50,
        pocket: PocketParameter {
            x: Vector3::new(0.0, 2.65e-3, 2.00e-3),
            r: 0.825e-3,
        },
        neck: NeckParameter {
            x: Vector3::new(0.0, 2.65e-3, 2.00e-3 - 0.30e-3),
            r: 1.20e-3 * 0.95,
            h: 2.45e-3 * 1.05,
            dh: 0.1519e-3 * 0.5,
            h_ratio: 0.355,
            r_ratio: 0.55,
        },
    };

    // let mut brg = Brg::new(&mesh, &param);
    let mut brg = Brg::new(&mesh, &param);

    let section0 = "edge";
    let section1 = "face";
    let edges = io::read_index_from_xml(index_path, section0).unwrap();
    let faces = io::read_index_from_xml(index_path, section1).unwrap();

    brg.set_edge_and_face(edges, faces);

    brg.linspace_all();

    brg.project_all();

    let iteration = 500;
    for _ in 0..iteration {
        brg.smooth_inner();
    }
    for _ in 0..iteration {
        brg.smooth_face();
    }
    for _ in 0..iteration {
        brg.smooth_inner();
    }
    for _ in 0..iteration {
        brg.smooth_face();
    }
    io::copy_vtk_and_replace_point(origin_path, write_path, &brg.get_points());

}




