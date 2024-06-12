pub mod vtk;
use morphing::vtk::brg::{Brg, CageParameter, };

use pyo3::prelude::*;
use CageParameter::*;

#[pymodule]
fn morphing(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<CageParameter>()?;
    Ok(())
}