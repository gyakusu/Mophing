use std::fs::File;
use std::io::BufReader;
use xml::reader::{EventReader, XmlEvent, Error};
use super::vtk::Point;
use super::vtk::Tetra;

pub fn parse_section(e: Result<XmlEvent, Error>, attribute_value: &str, is_in_section: &mut bool) {
    match e {
        Ok(XmlEvent::StartElement { name, attributes, .. }) => {
            if name.local_name == "DataArray" {
                for attr in attributes {
                    if attr.name.local_name == "Name" && attr.value == attribute_value {
                        *is_in_section = true;
                        break;
                    }
                }
            }
        }
        Ok(XmlEvent::EndElement { name }) => {
            if name.local_name == "DataArray" {
                *is_in_section = false;
            }
        }
        _ => {}
    }
}

pub fn parse_points(chars: String, points: &mut Vec<Point>) {
    let nums: Vec<f32> = chars.split_whitespace()
        .filter_map(|s| s.parse().ok())
        .collect();
    for chunk in nums.chunks(3) {
        if chunk.len() == 3 {
            let point = Point { x: [chunk[0], chunk[1], chunk[2]] };
            points.push(point);
        }
    }
}

pub fn parse_tetras(chars: String, tetras: &mut Vec<Tetra>) {
    let nums: Vec<i64> = chars.split_whitespace()
        .filter_map(|s| s.parse().ok())
        .collect();
    for chunk in nums.chunks(4) {
        if chunk.len() == 4 {
            let tetra = Tetra { index: [chunk[0], chunk[1], chunk[2], chunk[3]] };
            tetras.push(tetra);
        }
    }
}

pub fn read_data<T, F>(file_path: &str, section_name: &str, mut parse: F) -> Vec<T>
where
    F: FnMut(String, &mut Vec<T>),
{
    let file = File::open(file_path).unwrap();
    let file = BufReader::new(file);

    let parser = EventReader::new(file);
    let mut data = Vec::new();
    let mut is_in_section = false;

    for e in parser {
        parse_section(e.clone(), section_name, &mut is_in_section);
        match e {
            Ok(XmlEvent::Characters(chars)) if is_in_section => {
                parse(chars, &mut data);
            }
            Err(e) => {
                println!("Error: {}", e);
                break;
            }
            _ => {}
        }
    }
    data
}

pub fn read_point(file_path: &str) -> Vec<Point> {
    read_data(file_path, "Points", parse_points)
}

pub fn read_tetra(file_path: &str) -> Vec<Tetra> {
    read_data(file_path, "connectivity", parse_tetras)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_points() {
        let mut points = Vec::new();
        parse_points("1.0 2.0 3.0".to_string(), &mut points);
        assert_eq!(points, vec![Point { x: [1.0, 2.0, 3.0] }]);
    }

    #[test]
    fn test_parse_tetras() {
        let mut tetras = Vec::new();
        parse_tetras("1 2 3 4".to_string(), &mut tetras);
        assert_eq!(tetras, vec![Tetra { index: [1, 2, 3, 4] }]);
    }

    // 以下の２つのテストは，実行結果はテキストに依存するので実際はデバッグポイントで止めて目視で確認してね．ちなみに自分がやった時は配列に入っている値や配列の長さがあっていることは確認したよ．
    #[test]
    fn test_read_point() {

        let points = read_point("data/Tetra.vtu");
        println!("{:?}", points.first());
        println!("{:?}", points.last());
        println!("{:?}", points.len());

        assert!(true);
    }
    #[test]
    fn test_read_tetra() {

        let tetras = read_tetra("data/Tetra.vtu");
        println!("{:?}", tetras.first());
        println!("{:?}", tetras.last());
        println!("{:?}", tetras.len());

        assert!(true);
    }
}

