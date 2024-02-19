// main.rs
mod point;
mod cell;
mod face;
mod reader;

use point::Point;
use cell::Cell;
use face::Face;
use reader::Reader;

fn main() {
    let reader = Reader::new("path_to_your_file.vtk");
    let lines = reader.read_vtk_file().unwrap();

    let mut points = Vec::new();
    let mut cells = Vec::new();
    let mut faces = Vec::new();

    // main.rs
    for line in lines {

        // ここでは、行が "Points"、"Cells"、または "Faces" セクションの一部であるかどうかを判断します。
        // これは、VTKファイルの形式に依存します。
        if line.starts_with("Points") {
            
            // 点のデータを解析します。
            // ここでは、行からPoint構造体を作成し、pointsベクタに追加します。
            let point = parse_point(line);
            points.push(point);
        } 
        else if line.starts_with("Cells") {

            // セルのデータを解析します。
            // ここでは、行からCell構造体を作成し、cellsベクタに追加します。
            let cell = parse_cell(line);
            cells.push(cell);
        } 
        else if line.starts_with("Faces") {
        
            // 面のデータを解析します。
            // ここでは、行からFace構造体を作成し、facesベクタに追加します。
            let face = parse_face(line);
            faces.push(face);
        }
    }
}
