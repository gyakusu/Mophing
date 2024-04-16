use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::io::Read;
use std::io::Write;
use xml::reader::{EventReader, XmlEvent, Error};
use xml::writer;
use super::vtk::Face;
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
            let point = Point::new(chunk[0], chunk[1], chunk[2]);
            points.push(point);
        }
    }
}

pub fn parse_tetra(chars: String, tetras: &mut Vec<Tetra>) {
    let nums: Vec<usize> = chars.split_whitespace()
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
pub fn read_vtk(vtk_path: &str) -> Mesh {
    let points = read_point(vtk_path);
    let tetras = read_tetra(vtk_path);

    let mesh = Mesh::new(points, tetras);
    mesh
}
pub fn read_vtk_and_setting(vtk_path: &str, setting_path: &str) -> Mesh {
    let points = read_point(vtk_path);
    let tetras = read_tetra(vtk_path);

    let faces = read_face(setting_path, "face").unwrap();
    let surface_faces = read_face(setting_path, "surface_face").unwrap();
    let inner_index = read_index(setting_path, "inner_index").unwrap();
    let neighbor_map = read_map(setting_path, "neighbor_map").unwrap();
    let surface_map = read_map(setting_path, "surface_map").unwrap();

    Mesh::load(points, tetras, faces, surface_faces, inner_index, neighbor_map, surface_map)
}

pub fn write_setting(setting_path: &str, faces: Vec<Face>,  surface_faces: Vec<Face>, inner_index: Vec<usize>, neighbor_map: HashMap<usize, Vec<usize>>, surface_map: HashMap<usize, Vec<usize>>) -> std::io::Result<()> {
    let file = File::create(setting_path)?;
    let mut writer = writer::EmitterConfig::new()
        .perform_indent(true)
        .create_writer(file);

    let _ = writer.write(writer::XmlEvent::start_element("setting"));

    let _ = writer.write(writer::XmlEvent::start_element("face"));
    let _ = writer.write(writer::XmlEvent::characters("\n"));
    for face in faces {
        let face_str = face.as_i64().into_iter().map(|i| i.to_string()).collect::<Vec<String>>().join(" ");
        let _ = writer.write(writer::XmlEvent::characters(&face_str));
        let _ = writer.write(writer::XmlEvent::characters("\n"));
    }
    let _ = writer.write(writer::XmlEvent::end_element());

    let _ = writer.write(writer::XmlEvent::start_element("surface_face"));
    let _ = writer.write(writer::XmlEvent::characters("\n"));
    for face in surface_faces {
        let face_str = face.as_i64().into_iter().map(|i| i.to_string()).collect::<Vec<String>>().join(" ");
        let _ = writer.write(writer::XmlEvent::characters(&face_str));
        let _ = writer.write(writer::XmlEvent::characters("\n"));
    }
    let _ = writer.write(writer::XmlEvent::end_element());

    let _ = writer.write(writer::XmlEvent::start_element("inner_index"));
    let _ = writer.write(writer::XmlEvent::characters("\n"));
    for i in inner_index {
        let _ = writer.write(writer::XmlEvent::characters(&i.to_string()));
        let _ = writer.write(writer::XmlEvent::characters("\n"));
    }
    let _ = writer.write(writer::XmlEvent::end_element());

    let _ = writer.write(writer::XmlEvent::start_element("neighbor_map"));
    let _ = writer.write(writer::XmlEvent::characters("\n"));
    for (key, values) in neighbor_map {
        let values_str = values.into_iter().map(|v| v.to_string()).collect::<Vec<String>>().join(" ");
        let line = format!("{} {}", key, values_str);
        let _ = writer.write(writer::XmlEvent::characters(&line));
        let _ = writer.write(writer::XmlEvent::characters("\n"));
    }
    let _ = writer.write(writer::XmlEvent::end_element());

    let _ = writer.write(writer::XmlEvent::start_element("surface_map"));
    let _ = writer.write(writer::XmlEvent::characters("\n"));
    for (key, values) in surface_map {
        let values_str = values.into_iter().map(|v| v.to_string()).collect::<Vec<String>>().join(" ");
        let line = format!("{} {}", key, values_str);
        let _ = writer.write(writer::XmlEvent::characters(&line));
        let _ = writer.write(writer::XmlEvent::characters("\n"));
    }
    let _ = writer.write(writer::XmlEvent::end_element());

    let _ = writer.write(writer::XmlEvent::end_element());

    Ok(())
}

pub fn read_face(setting_path: &str, target: &str) -> std::io::Result<Vec<Face>> {
    let file = File::open(setting_path)?;
    let file = BufReader::new(file);

    let parser = EventReader::new(file);
    let mut faces = Vec::new();
    let mut current_face = Vec::new();

    for e in parser {
        match e {
            Ok(XmlEvent::StartElement { name, .. }) if name.local_name == target => {
                current_face.clear();
            }
            Ok(XmlEvent::Characters(s)) => {
                let numbers: Vec<usize> = s.split_whitespace().map(|s| s.parse().unwrap()).collect();
                current_face.extend(numbers);
            }
            Ok(XmlEvent::EndElement { name }) if name.local_name == target => {
                while current_face.len() >= 3 {
                    let face = Face::new([current_face[0], current_face[1], current_face[2]]).unwrap();
                    faces.push(face);
                    current_face.drain(0..3);
                }
            }
            Err(e) => {
                println!("Error: {}", e);
                break;
            }
            _ => {}
        }
    }
    Ok(faces)
}

pub fn read_index(setting_path: &str, target: &str) -> std::io::Result<Vec<usize>> {
    let file = File::open(setting_path)?;
    let file = BufReader::new(file);

    let parser = EventReader::new(file);
    let mut index = Vec::new();
    let mut is_target_section = false;

    for e in parser {
        match e {
            Ok(XmlEvent::StartElement { name, .. }) => {
                is_target_section = name.local_name == target;
            }
            Ok(XmlEvent::Characters(s)) if is_target_section => {
                let numbers: Vec<usize> = s.split_whitespace()
                                         .filter_map(|s| s.parse().ok())
                                         .collect();
                index.extend(numbers);
            }
            Ok(XmlEvent::EndElement { name }) if name.local_name == target => {
                is_target_section = false;
            }
            Err(e) => {
                println!("Error: {}", e);
                break;
            }
            _ => {}
        }
    }
    Ok(index)
}

pub fn read_map(setting_path: &str, target: &str) -> std::io::Result<HashMap<usize, Vec<usize>>> {
    let file = File::open(setting_path)?;
    let file = BufReader::new(file);

    let parser = EventReader::new(file);
    let mut neighbor_map = HashMap::new();
    let mut is_target_section = false;

    for e in parser {
        match e {
            Ok(XmlEvent::StartElement { name, .. }) => {
                is_target_section = name.local_name == target;
            }
            Ok(XmlEvent::Characters(s)) if is_target_section => {
                for line in s.lines() {
                    let numbers: Vec<usize> = line.split_whitespace()
                                                .filter_map(|s| s.parse().ok())
                                                .collect();
                    if let Some(key) = numbers.get(0) {
                        let values = numbers[1..].to_vec();
                        neighbor_map.insert(*key, values);
                    }
                }
            }
            Ok(XmlEvent::EndElement { name }) if name.local_name == target => {
                is_target_section = false;
            }
            Err(e) => {
                println!("Error: {}", e);
                break;
            }
            _ => {}
        }
    }
    Ok(neighbor_map)
}
pub fn read_index_from_xml(file_path: &str, section: &str) -> std::io::Result<HashMap<String, Vec<usize>>> {

    let file = File::open(file_path)?;
    let file = BufReader::new(file);

    let parser = EventReader::new(file);
    let mut edges = HashMap::new();
    let mut current_edge = String::new();
    let mut in_edge_section = false;
    
    for e in parser {
        match e {
            Ok(XmlEvent::StartElement { name, .. }) => {
                if name.local_name == section {
                    in_edge_section = true;
                } else if in_edge_section {
                    current_edge = name.local_name;
                }
            }
            Ok(XmlEvent::EndElement { name, .. }) => {
                if name.local_name == section {
                    in_edge_section = false;
                }
            }
            Ok(XmlEvent::Characters(s)) => {
                if in_edge_section {
                    let numbers: Vec<usize> = s.split_whitespace()
                                            .filter_map(|s| s.parse().ok())
                                            .collect();
                    edges.insert(current_edge.clone(), numbers);
                }
            }
            Err(e) => {
                println!("Error: {}", e);
                break;
            }
            _ => {}
        }
    }
    Ok(edges)
}
pub fn copy_vtk_and_replace_point(old_file_path: &str, new_file_path: &str, points: &Vec<Point>) {
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
        let x = point.as_vec();
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
        assert_eq!(points, vec![Point::new(1.0, 2.0, 3.0)] );
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
    fn test_read_vtk() {
        let mesh = read_vtk("data/Tetra.vtu");
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
            Point::new(1.0, 2.0, 3.0),
            Point::new(4.0, 5.0, 6.0),
        ];
        copy_vtk_and_replace_point("data/Tetra.vtu", "data/Tetra_copy.vtu", &mut points);
        assert!(true);
    }
    #[test]
    fn test_read_index_from_xml() {
        let file_path = "data/face_and_edge_index.xml";
        let section0 = "edge";
        let section1 = "face";
        let edges = read_index_from_xml(file_path, section0).unwrap();
        let faces = read_index_from_xml(file_path, section1).unwrap();
        // 目視で確認してください．
        assert!(edges.len() > 0);
        assert!(faces.len() > 0);
    }
    #[test]
    #[ignore]
    fn test_write_and_read_setting() {
        let mesh = read_vtk("data/Tetra.vtu");

        // mesh.save();
        let(points, tetras, faces, surface_faces, inner_index, neighbor_map, surface_map) = mesh.save();
        
        _ = write_setting("data/Tetra_setting.xml", faces.clone(), surface_faces.clone(), inner_index.clone(), neighbor_map.clone(), surface_map.clone());

        let faces0 = read_face("data/Tetra_setting.xml", "face").unwrap();
        let surface_faces0 = read_face("data/Tetra_setting.xml", "surface_face").unwrap();
        let inner_index0 = read_index("data/Tetra_setting.xml", "inner_index").unwrap();
        let neighbor_map0 = read_map("data/Tetra_setting.xml", "neighbor_map").unwrap();
        let surface_map0 = read_map("data/Tetra_setting.xml", "surface_map").unwrap();

        let _ = Mesh::load(points, tetras, faces0.clone(), surface_faces0.clone(), inner_index0.clone(), neighbor_map0.clone(), surface_map0.clone());

        assert_eq!(faces.clone().first(), faces0.clone().first());
        assert_eq!(faces.clone().last(), faces0.clone().last());

        assert_eq!(surface_faces.clone().first(), surface_faces0.clone().first());
        assert_eq!(surface_faces.clone().last(), surface_faces0.clone().last());

        assert_eq!(inner_index.clone().first(), inner_index0.clone().first());
        assert_eq!(inner_index.clone().last(), inner_index0.clone().last());

        assert_eq!(neighbor_map.clone().get(&1).unwrap().first(), neighbor_map0.clone().get(&1).unwrap().first());
        assert_eq!(neighbor_map.clone().get(&1).unwrap().last(), neighbor_map0.clone().get(&1).unwrap().last());

        assert_eq!(surface_map.clone().get(&10).unwrap().first(), surface_map0.clone().get(&10).unwrap().first());
        assert_eq!(surface_map.clone().get(&10).unwrap().last(), surface_map0.clone().get(&10).unwrap().last());
    }
    #[test]
    #[ignore]
    fn test_read_vtk_and_setting() {
        let mesh0 = read_vtk_and_setting("data/Tetra.vtu", "data/Tetra_setting.xml");
        let mesh1 = read_vtk("data/Tetra.vtu");

        assert_eq!(mesh0.points.last(), mesh1.points.last());
        assert_eq!(mesh0.inner_index.last(), mesh1.inner_index.last());
        assert_eq!(mesh0.faces.last(), mesh1.faces.last());
        assert_eq!(mesh0.tetras.last(), mesh1.tetras.last());
        assert_eq!(mesh0.neighbor_map.get(&1).unwrap().last(), mesh1.neighbor_map.get(&1).unwrap().last());
        assert_eq!(mesh0.surface_map.get(&10).unwrap().last(), mesh1.surface_map.get(&10).unwrap().last());

        assert!(true);
    }
}

