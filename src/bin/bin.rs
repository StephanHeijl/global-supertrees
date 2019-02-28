#[macro_use(s)]
extern crate ndarray;
extern crate rayon;
extern crate petgraph;
extern crate regex;
extern crate serde;
extern crate bincode;
extern crate global_supertrees;


use std::fs::File;
use std::io::prelude::*;
use std::env;
use rayon::prelude::*;
use global_supertrees::tree_distance_matrix::TreeDistanceMatrix;
use std::collections::BTreeMap;
use global_supertrees::tree_merging::merge_trees;



fn main() {

    let mut mapping : BTreeMap<[u8; 10], u32> = BTreeMap::new();
    let mut trees : Vec<String> = Vec::new();

    for arg in env::args() {
        if arg.ends_with(".txt") {
            let cache_mapping_path = global_supertrees::utils::get_cache_mapping_name(&arg);
            match global_supertrees::utils::load_cache_mapping(&cache_mapping_path) {
                Ok(m) => {
                    mapping = m;
                    println!("Used the cached the mapping file here: {:?} ", cache_mapping_path);
                },
                Err(_e) => {
                    mapping = global_supertrees::uniprot::load_mapping_file(&arg);
                    global_supertrees::utils::cache_mapping(&mapping, &cache_mapping_path).expect(
                        "Could not cache mapping file."
                    );
                    println!("Cached the mapping file here: {:?} ", cache_mapping_path);
                }
            }

        } else if arg.ends_with(".tree") || arg.ends_with(".ftree") {
            trees.push(arg);
        }
    }

    let mut distance_matrices : Vec<TreeDistanceMatrix> = trees.par_iter().map(
        |t| global_supertrees::utils::convert_file_to_distance_matrix(t.to_string())
    ).collect();

    if mapping.len() > 0 {
        distance_matrices = distance_matrices.par_iter().map(
            |dm| dm.merge_organisms(&mapping)
        ).collect();
    }

    let trees : Vec<global_supertrees::graph_tree::Tree> = distance_matrices.par_iter().map(
        | dm | dm.neighbour_joining()
    ).collect();

    for (t, tree) in trees.iter().enumerate() {
        let mut file = File::create(format!("sub_{}.tree", t)).expect("IO Error while creating file.");
        file.write_all(tree.to_newick().as_bytes()).expect("IO Error while writing merged tree.");
    }

    let merged_tree = merge_trees(trees);

    let mut file = File::create("merged.tree").expect("IO Error while creating file.");
    file.write_all(merged_tree.to_newick().as_bytes()).expect("IO Error while writing merged tree.");


}
