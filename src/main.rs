use nalgebra as na;
use na::Vector3;
use std::f32::consts::PI;

use mophing::vtk::io;
use mophing::vtk::brg::{Brg, CageParameter, PocketParameter, NeckParameter};
use mophing::vtk::vtk::check_laplacian_smoothing_quality;

// example: `cargo run --release "data/Tetra.vtu,data/face_and_edge_index.xml,data/Tetra_linspace.vtu" "2.345e-3,2.850e-3,0.93e-3,2.10e-3,0.10e-3,0.825e-3,1.70e-3,1.20e-3,2.45e-3,0.152e-3" "100"`
fn main() {

    let args: Vec<String> = std::env::args().collect();
    let parsed_args: Vec<Vec<String>> = args[1..]
        .iter()
        .map(|arg| {
            arg.split(',')
                .map(|s| s.to_string())
                .collect()
        })
        .collect();

    let get_str = |i: usize| parsed_args[0][i].as_str();
    let get_f32 = |i: usize| parsed_args[1][i].parse::<f32>().unwrap();
    let get_usize = |i: usize| parsed_args[2][i].parse::<usize>().unwrap();

    println!("{:?}", parsed_args);

    let origin_path = get_str(0);
    let index_path = get_str(1);
    let write_path = get_str(2);

    let mesh = io::read_vtk(origin_path);

    let param = CageParameter {
        axis: Vector3::new(0.0, 0.0, 1.0),
        theta0: -PI / 6.0,
        theta1:  PI / 6.0,
        r0: get_f32(0),
        r1: get_f32(1),
        h0: get_f32(2),
        h1: get_f32(3),
        bevel: get_f32(4),
        pocket: PocketParameter {
            x: Vector3::new(0.0, 2.65e-3, 2.00e-3),
            r: get_f32(5),
        },
        neck: NeckParameter {
            x: Vector3::new(0.0, 2.65e-3, get_f32(6)),
            r: get_f32(7),
            h: get_f32(8),
            dh: get_f32(9),
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

    let iteration = get_usize(0);
    
    for i in 0..iteration {
        let old_points = brg.get_points();
        brg.smooth_face();
        let new_points = brg.get_points();

        if i != 0 && (i & (i - 1)) == 0 {
            let quality = check_laplacian_smoothing_quality(old_points, new_points);
            print!("iteration: {:?}, ", i);
            println!("quality: {:?}", quality);
            if quality < 1e-7 {
                break;
            }
        }
    }
    for i in 0..iteration {
        let old_points = brg.get_points();
        brg.smooth_inner();
        let new_points = brg.get_points();

        if i != 0 && (i & (i - 1)) == 0 {
            let quality = check_laplacian_smoothing_quality(old_points, new_points);
            print!("iteration: {:?}, ", i);
            println!("quality: {:?}", quality);
            if quality < 1e-7 {
                break;
            }
        }
    }
    io::copy_vtk_and_replace_point(origin_path, write_path, &brg.get_points());

}




