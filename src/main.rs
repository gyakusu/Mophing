use morphing::vtk::brg::{Brg, CageParameter, };
use morphing::vtk::mesh::check_smoothing_quality;

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

    let cage = CageParameter::new(arg_f32(0), arg_f32(1), arg_f32(2), arg_f32(3), arg_f32(4), arg_f32(5), arg_f32(6), arg_f32(7), arg_f32(8), arg_f32(9));
    let mut brg = Brg::from_file(origin_path, index_path, &cage);

    brg.linspace_all();

    let iterations = vec![arg_usize(0), arg_usize(1), arg_usize(2), arg_usize(3)];
    
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
    smooth_and_check_quality(&mut brg, iterations[0], |brg| {
        brg.smooth_face()
    });
    smooth_and_check_quality(&mut brg, iterations[1], |brg| {
        brg.smooth_ball();
        brg.smooth_inner();
    });

    let scale = brg.std();
    println!("scale: {:?}", scale);

    println!("flip negative volume start");
    for _ in 0..iterations[2] {
        brg.scale(1.0/scale);
        let success = brg.flip_negative_volume(0.1);
        brg.scale(scale);
        brg.normalize_center();
        if success {
            break;
        }
    }

    println!("flip negative volume only sphire start");
    let mut all_success = false;
    for _ in 0..iterations[3] {
        brg.scale(1.0/scale);
        let success = brg.flip_negative_volume_only_sphire(2.0);
        brg.scale(scale);
        brg.normalize_center();
        if success {
            all_success = true;
            break;
        }
    }
    if all_success {
        println!("I finished the all process.");
        brg.write_vtk_from_base(origin_path, write_path);
    }
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




