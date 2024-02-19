use mophing::reader::Reader;
use std::fs::File;
use std::io::Write;
use std::path::Path;
use std::io::BufReader;
use xml::reader::{EventReader, XmlEvent};

#[test]
fn test_read_vtk_file() {
    // テスト用のVTKファイルを作成します。
    let test_file_path = "test.vtk";
    let mut file = File::create(&test_file_path).unwrap();
    writeln!(file, "line1\nline2\nline3").unwrap();

    // Readerを作成し、テスト用のファイルを読み込みます。
    let reader = Reader::new(test_file_path);
    let lines = reader.read_vtk_file().unwrap();

    // 読み込んだ行が期待通りであることを確認します。
    assert_eq!(lines.len(), 3);
    assert_eq!(lines[0], "line1");
    assert_eq!(lines[1], "line2");
    assert_eq!(lines[2], "line3");

    // テスト用のファイルを削除します。
    std::fs::remove_file(Path::new(test_file_path)).unwrap();
}

#[test]
fn test_read_vtk_file2() {
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
