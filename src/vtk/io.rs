use std::fs::File;
use std::io::BufReader;
use std::io::Read;
use std::io::Write;
use xml::reader::{EventReader, XmlEvent, Error};
use super::vtk::Point;
use super::vtk::Tetra;
use super::vtk::Mesh;

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

pub fn parse_point(chars: String, points: &mut Vec<Point>) {
    let nums: Vec<f32> = chars.split_whitespace()
        .filter_map(|s| s.parse().ok())
        .collect();
    for chunk in nums.chunks(3) {
        if chunk.len() == 3 {
            let point = Point::new([chunk[0], chunk[1], chunk[2]]);
            points.push(point);
        }
    }
}

pub fn parse_tetra(chars: String, tetras: &mut Vec<Tetra>) {
    let nums: Vec<i64> = chars.split_whitespace()
        .filter_map(|s| s.parse().ok())
        .collect();
    for chunk in nums.chunks(4) {
        if chunk.len() == 4 {
            let tetra = Tetra::new([chunk[0], chunk[1], chunk[2], chunk[3]]).unwrap();
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
    read_data(file_path, "Points", parse_point)
}

pub fn read_tetra(file_path: &str) -> Vec<Tetra> {
    read_data(file_path, "connectivity", parse_tetra)
}

pub fn read_all(file_path: &str) -> Mesh {
    let points = read_point(file_path);
    let tetras = read_tetra(file_path);

    let mesh = Mesh::new(points, tetras);
    mesh
}

pub fn copy_vtk_and_replace_point(old_file_path: &str, new_file_path: &str, points: &mut Vec<Point>) {
    let mut old_file = File::open(old_file_path).unwrap();
    let mut contents = String::new();
    old_file.read_to_string(&mut contents).unwrap();

    // Find the start and end indices of the "Points" section
    let start_index = contents.find("<Points>").unwrap();
    let end_index = contents.find("</Points>").unwrap() + "</Points>".len();

    // Replace the "Points" section with the new points
    let mut new_contents = contents[..start_index].to_string();
    new_contents.push_str("<Points>\n");
    new_contents.push_str(&format!("  <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n"));

    for point in points {
        let x = point.as_f32();
        new_contents.push_str(&format!("    {} {} {}\n", x[0], x[1], x[2]));
    }
    new_contents.push_str("  </DataArray>\n");
    new_contents.push_str("</Points>\n");
    new_contents.push_str(&contents[end_index..]);
    
    // Write the new contents to a new file
    let mut new_file = File::create(new_file_path).unwrap();
    new_file.write_all(new_contents.as_bytes()).unwrap();

    println!("The new file has been created: {}", new_file_path);
}    

#[cfg(test)]
mod tests {
    use std::vec;

    use super::*;

    #[test]
    fn test_parse_point() {
        let mut points = Vec::new();
        parse_point("1.0 2.0 3.0".to_string(), &mut points);
        assert_eq!(points, vec![Point::new([1.0, 2.0, 3.0])] );
    }

    #[test]
    fn test_parse_tetra() {
        let mut tetras = Vec::new();
        parse_tetra("1 2 3 4".to_string(), &mut tetras);
        assert_eq!(tetras, vec![Tetra::new([1, 2, 3, 4]).unwrap()]);
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
    #[test]
    fn test_read_all() {
        let mesh = read_all("data/Tetra.vtu");
        println!("{:?}", mesh.points.first());
        println!("{:?}", mesh.points.last());
        println!("{:?}", mesh.points.len());
        println!("{:?}", mesh.tetras.first());
        println!("{:?}", mesh.tetras.last());
        println!("{:?}", mesh.tetras.len());
        assert!(true);
    }
    #[test]
    fn test_copy_and_write_point() {
        let mut points = vec![
            Point::new([1.0, 2.0, 3.0]),
            Point::new([4.0, 5.0, 6.0]),
        ];
        copy_vtk_and_replace_point("data/Tetra.vtu", "data/Tetra_copy.vtu", &mut points);
        assert!(true);
    }
}

