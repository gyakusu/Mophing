pub mod vtk;
use vtk::brg::{CageParameter, Brg};

use pyo3::prelude::*;

#[pymodule]
fn morphing(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<CageParameter>()?;
    m.add_class::<Brg>()?;
    Ok(())
}


