// use std::collections::HashMap;

// // 座標を表す構造体
// #[derive(Debug, Clone)]
// struct Point {
//     x: f64,
//     y: f64,
//     z: f64,
// }

// // 結合情報を表す構造体
// #[derive(Debug, Clone)]
// struct Connection {
//     point1: Point,
//     point2: Point,
// }

// // テトラメッシュを表す構造体
// #[derive(Debug)]
// struct TetraMesh {
//     points: HashMap<String, Point>,
//     connections: HashMap<String, Connection>,
// }

// fn main() {
//     // ポイントと接続を作成
//     let point1 = Point { x: 0.0, y: 0.0, z: 0.0 };
//     let point2 = Point { x: 1.0, y: 0.0, z: 0.0 };
//     let connection = Connection { point1: point1.clone(), point2: point2.clone() };

//     // ポイントと接続をハッシュマップに追加
//     let mut points = HashMap::new();
//     points.insert("point1".to_string(), point1);
//     points.insert("point2".to_string(), point2);

//     let mut connections = HashMap::new();
//     connections.insert("connection1".to_string(), connection);

//     // テトラメッシュを作成
//     let tetra_mesh = TetraMesh { points, connections };

//     // テトラメッシュの内容を表示
//     println!("{:#?}", tetra_mesh);
// }

use std::fs::File;
use std::io::BufReader;
use xml::reader::{EventReader, XmlEvent};

fn main() {
    let file = File::open("data/Trim_Cells9291_0_0.vtu").unwrap();
    let file = BufReader::new(file);

    let parser = EventReader::new(file);
    for e in parser {
        match e {
            Ok(XmlEvent::StartElement { name, attributes, .. }) => {
                if name.local_name == "DataArray" {
                    println!("Attributes of DataArray:");
                    for attr in attributes {
                        println!("{} = {}", attr.name.local_name, attr.value);
                    }
                }
            }
            Err(e) => {
                println!("Error: {}", e);
                break;
            }
            _ => {}
        }
    }
}
