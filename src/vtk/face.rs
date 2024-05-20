extern crate nalgebra as na;

#[derive(Debug, Eq, Hash, PartialEq, Clone, Copy)]
pub struct Face {
    pub index: [usize; 3],
}

impl Face {
    pub fn new(index: [usize; 3]) -> Result<Self, &'static str> {
        if index.len() != 3 {
            return Err("index must have a length of 3");
        }        
        let mut sorted_index = index;
        sorted_index.sort();
        Ok(Self {
            index: sorted_index,
        })
    }
    pub fn as_i64(&self) -> [usize; 3] {
        self.index
    }
}

