// reader.rs
use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;

pub struct Reader {
    pub filename: String,
}

impl Reader {
    pub fn new(filename: &str) -> Reader {
        Reader {
            filename: filename.to_string(),
        }
    }

    pub fn read_vtk_file(&self) -> io::Result<Vec<String>> {
        let path = Path::new(&self.filename);
        let file = File::open(&path)?;
        let reader = io::BufReader::new(file);

        let mut lines = Vec::new();

        for line in reader.lines() {
            let line = line?;
            lines.push(line);
        }

        Ok(lines)
    }
}
