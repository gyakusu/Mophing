use nalgebra as na;
use na::Vector3;
use std::f32::consts::PI;

use mophing::vtk::io;
use mophing::vtk::brg::{Brg, CageParameter, PocketParameter, NeckParameter};
use mophing::vtk::mesh::check_smoothing_quality;

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

    let arg_str = |i: usize| parsed_args[0][i].as_str();
    let arg_f32 = |i: usize| parsed_args[1][i].parse::<f32>().unwrap();
    let arg_usize = |i: usize| parsed_args[2][i].parse::<usize>().unwrap();

    println!("{:?}", parsed_args);

    let origin_path = arg_str(0);
    let index_path = arg_str(1);
    let write_path = arg_str(2);

    let mesh = io::read_vtk(origin_path);

    println!("points: {:?}", mesh.points.len());
    println!("tetras: {:?}", mesh.tetras.len());
    println!("inner_index: {:?}", mesh.inner_index.len());
    println!("inverse_map: {:?}", mesh.inverse_map.len());


    let param = CageParameter {
        axis: Vector3::new(0.0, 0.0, 1.0),
        theta0: -PI / 6.0,
        theta1:  PI / 6.0,
        r0: arg_f32(0),
        r1: arg_f32(1),
        h0: arg_f32(2),
        h1: arg_f32(3),
        bevel: arg_f32(4),
        pocket: PocketParameter {
            x: Vector3::new(0.0, 2.65e-3, 2.00e-3),
            r: arg_f32(5),
        },
        neck: NeckParameter {
            x: Vector3::new(0.0, 2.65e-3, arg_f32(6)),
            r: arg_f32(7),
            h: arg_f32(8),
            dh: arg_f32(9),
            h_ratio: 0.355,
            r_ratio: 0.55,
        },
    };

    let mut brg = Brg::new(&mesh, &param);

    let section0 = "edge";
    let section1 = "face";
    let edges = io::read_edge_from_xml(index_path, section0).unwrap();
    let faces = io::read_index_from_xml(index_path, section1).unwrap();

    brg.set_edge_and_face(&edges, &faces);

    println!("set edge and face done");
    println!("edge: {:?}", edges.len());

    brg.linspace_all();
    brg.project_all();

    let iteration0 = arg_usize(0);
    let iteration1 = arg_usize(1);
    let iteration2 = arg_usize(2);
    let iteration3 = arg_usize(3);
    
    fn smooth_and_check_quality(brg: &mut Brg, iteration: usize, smooth_fn: impl Fn(&mut Brg)) {
        for i in 0..iteration {
            let old_points = brg.get_points();
            smooth_fn(brg);
    
            if i != 0 && (i & (i - 1)) == 0 {
                let new_points = brg.get_points();
                let quality = check_smoothing_quality(&old_points, &new_points);
                print!("iteration: {:?}, ", i);
                println!("quality: {:?}", quality);
                if quality < 1e-9 {
                    break;
                }
            }
        }
    }
    smooth_and_check_quality(&mut brg, iteration0, |brg| {
        brg.smooth_face()
    });
    smooth_and_check_quality(&mut brg, iteration1, |brg| {
        brg.smooth_ball();
        brg.smooth_inner();
    });

    let scale = brg.std();
    println!("scale: {:?}", scale);

    println!("flip negative volume start");
    for _ in 0..iteration2 {
        brg.scale(1.0/scale);
        brg.flip_negative_volume(0.1);

        brg.scale(scale);
        brg.normalize_center();
    }
    println!("flip negative volume only sphire start");
    for _ in 0..iteration3 {
        brg.scale(1.0/scale);
        brg.flip_negative_volume_only_sphire(2.0);

        brg.scale(scale);
        brg.normalize_center();
    }
    io::copy_vtk_and_replace_point(origin_path, write_path, &brg.get_points());
}


    // brg.linspace_ball();
    // println!("linspace ball done");
    // for i in 0..iteration {
    //     let old_points = brg.get_points();
    //     brg.smooth_inner();

    //     if i != 0 && (i & (i - 1)) == 0 {
    //         let new_points = brg.get_points();
    //         let quality = check_smoothing_quality(old_points, new_points);
    //         print!("iteration: {:?}, ", i);
    //         println!("quality: {:?}", quality);
    //         if quality < 1e-14 {
    //             break;
    //         }
    //     }
    // }
    // let ball_map = io::read_map(index_path, "ball_map").unwrap();
    // brg.set_ball_map(ball_map);




