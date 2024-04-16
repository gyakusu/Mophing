use nalgebra as na;
use na::Vector3;
use std::f32::consts::PI;

use mophing::vtk::io;
use mophing::vtk::brg::Brg;
use mophing::vtk::brg::{CageParameter, PocketParameter, NeckParameter};

// example: `cargo run --release -- data/Tetra_linspace.vtu 2.345e-3 2.850e-3 0.93e-3 2.10e-3 0.10e-3 0.825e-3 1.70e-3 1.20e-3 2.45e-3 0.152e-3 100`
fn main() {

    let args: Vec<String> = std::env::args().collect();
    let get_arg = |i: usize| args[i].parse::<f32>().unwrap();

    println!("Argument:");
    println!("{}", args.join(" "));

    let origin_path = "data/Tetra.vtu";
    let index_path = "data/face_and_edge_index.xml";
    let write_path = args[1].as_str();

    let mesh = io::read_vtk(origin_path);

    let param = CageParameter {
        axis: Vector3::new(0.0, 0.0, 1.0),
        theta0: -PI / 6.0,
        theta1:  PI / 6.0,
        r0: get_arg(2),
        r1: get_arg(3),
        h0: get_arg(4),
        h1: get_arg(5),
        bevel: get_arg(6),
        pocket: PocketParameter {
            x: Vector3::new(0.0, 2.65e-3, 2.00e-3),
            r: get_arg(7),
        },
        neck: NeckParameter {
            x: Vector3::new(0.0, 2.65e-3, get_arg(8)),
            r: get_arg(9),
            h: get_arg(10),
            dh: get_arg(11),
            h_ratio: 0.355,
            r_ratio: 0.55,
        },
    };

    let mut brg = Brg::new(&mesh, &param);

    let section0 = "edge";
    let section1 = "face";
    let edges = io::read_index_from_xml(index_path, section0).unwrap();
    let faces = io::read_index_from_xml(index_path, section1).unwrap();

    brg.set_edge_and_face(edges, faces);

    brg.linspace_all();

    brg.project_all();

    let iteration = args[12].parse::<usize>().unwrap();
    for _ in 0..iteration {
        brg.smooth_face();
    }
    for _ in 0..iteration {
        brg.smooth_inner();
    }
    io::copy_vtk_and_replace_point(origin_path, write_path, &brg.get_points());

}




