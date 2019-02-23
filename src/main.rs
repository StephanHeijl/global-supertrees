#[macro_use(s)]
extern crate ndarray;
extern crate rayon;
extern crate petgraph;
extern crate regex;
extern crate serde;
extern crate bincode;

use std::env;
use std::time;

use rayon::prelude::*;
use tree_distance_matrix::TreeDistanceMatrix;
use std::collections::BTreeMap;

mod test_tree_parsing;
mod tests;
pub mod graph_tree;
pub mod tree_distance_matrix;
pub mod tree_merging;
pub mod utils;
pub mod uniprot;


pub fn convert_file_to_distance_matrix(fname : String) -> TreeDistanceMatrix {
    let now = time::Instant::now();
    let tree_file = utils::load_tree_file(String::from(fname));
    println!(
        "Loaded tree in {}.{} seconds.",
        now.elapsed().as_secs(),
        now.elapsed().subsec_millis()
    );

    let now = time::Instant::now();
    let parsed_tree = graph_tree::Tree::parse(tree_file);
    println!(
        "Parsed tree in {}.{} seconds.",
        now.elapsed().as_secs(),
        now.elapsed().subsec_millis()
    );

    let now = time::Instant::now();
    let _children = parsed_tree.traverse_children();
    println!(
        "Built tree traversal map in {}.{} seconds.",
        now.elapsed().as_secs(),
        now.elapsed().subsec_millis()
    );

    let now = time::Instant::now();
    let distance_matrix = parsed_tree.to_distance_matrix();
    println!(
        "Built distance_matrix in {}.{} seconds.",
        now.elapsed().as_secs(),
        now.elapsed().subsec_millis()
    );

    return distance_matrix;
}


fn main() {

    let mut mapping : BTreeMap<[u8; 10], u32> = BTreeMap::new();

    let mut trees : Vec<String> = Vec::new();

    for arg in env::args() {
        if arg.ends_with(".txt") {
            let cache_mapping_path = utils::get_cache_mapping_name(&arg);
            match utils::load_cache_mapping(&cache_mapping_path) {
                Ok(m) => {
                    mapping = m;
                    println!("Used the cached the mapping file here: {:?} ", cache_mapping_path);
                },
                Err(_e) => {
                    mapping = uniprot::load_mapping_file(&arg);
                    utils::cache_mapping(&mapping, &cache_mapping_path).expect("Could not cache mapping file.");
                    println!("Cached the mapping file here: {:?} ", cache_mapping_path);
                }
            }

        } else if arg.ends_with(".tree") || arg.ends_with(".ftree") {
            trees.push(arg);
        }
    }

    let mut distance_matrices : Vec<TreeDistanceMatrix> = trees.par_iter().map(
        |t| convert_file_to_distance_matrix(t.to_string())
    ).collect();

    if mapping.len() > 0 {
        distance_matrices = distance_matrices.par_iter().map(
            |dm| dm.merge_organisms(&mapping)
        ).collect();
    }

    let trees : Vec<graph_tree::Tree> = distance_matrices.par_iter().map(| dm | dm.neighbour_joining()).collect();

    tree_merging::merge_trees(trees);
}
