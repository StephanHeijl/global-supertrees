extern crate pyo3;
#[macro_use(s)]
extern crate ndarray;
extern crate numpy;
extern crate rayon;
extern crate petgraph;
extern crate regex;
extern crate serde;
extern crate bincode;

pub mod test_tree_parsing;
pub mod tests;

pub mod graph_tree;
pub mod tree_distance_matrix;
pub mod tree_merging;
pub mod utils;
pub mod uniprot;

use pyo3::prelude::*;
use ndarray::Array2;
use utils::convert_file_to_distance_matrix;

use numpy::{IntoPyArray, PyArray2};


#[pymodinit]
fn rust_ext(_py: Python, m: &PyModule) -> PyResult<()> {
    fn convert_file_to_dm(fname : String) -> Array2<f32> {
        let dm = convert_file_to_distance_matrix(fname);
        return dm.distance_matrix;
    }

    // wrapper of `axpy`
    #[pyfn(m, "convert_file_to_distance_matrix")]
    fn convert_file_to_distance_matrix_py(py: Python, fname: String) -> Py<PyArray2<f32>> {
        convert_file_to_dm(fname).into_pyarray(py).to_owned()
    }


    Ok(())
}