// vtk_parser.rs
use std::fs::File;
use std::io::BufReader;
use xml::reader::{EventReader, XmlEvent};
use super::point::Point;

pub fn parse_array(file_path: &str, section_keyword: &str) -> Vec<Point> {
    let file = File::open(file_path).unwrap();
    let file = BufReader::new(file);

    let parser = EventReader::new(file);
    let mut points = Vec::new();

    for e in parser {
        match e {
            Ok(XmlEvent::StartElement { name, attributes, .. }) => {
                if name.local_name == section_keyword {
                    // ここでPoint配列を作成し、それをpointsに追加します。
                    // Pointの定義とその作成方法は、VTKファイルの形式と内容によります。
                    // 以下は一例です。
                    let parse_attribute = |name: &str| {
                        attributes.iter()
                            .find(|a| a.name.local_name == name)
                            .and_then(|a| a.value.parse::<f32>().ok())
                            .unwrap_or(0.0)
                    };
                    
                    let point = Point {
                        x: parse_attribute("x"),
                        y: parse_attribute("y"),
                        z: parse_attribute("z"),
                    };                    
                    points.push(point);
                }
            }
            Err(e) => {
                println!("Error: {}", e);
                break;
            }
            _ => {}
        }
    }

    points
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_my_function() {
        let points = parse_array("data/Trim_Cells9291_0_0.vtu", "DataArray");
        println!("{:?}", points);
    }
}

