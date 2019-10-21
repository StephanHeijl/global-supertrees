#[macro_use(s)]
extern crate ndarray;
extern crate rayon;
extern crate petgraph;
extern crate regex;
extern crate serde;
extern crate bincode;
extern crate rand;

use std::f32;
use std::fs::File;
use std::io::prelude::*;
use std::env;
use rayon::prelude::*;
use rand::prelude::*;
use std::collections::BTreeMap;
use crate::tree_merging::mean_merge_distance_matrices;
use crate::graph_tree::Tree;
use crate::tree_distance_matrix::TreeDistanceMatrix;


fn main() {
    let mut mapping : BTreeMap<[u8; 10], u32> = BTreeMap::new();
    let mut trees : Vec<String> = Vec::new();

    let mut command : String = String::from("merge");

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

        } else if arg.ends_with(".tree") || arg.ends_with(".ftree") || arg.ends_with(".tre") {
            trees.push(arg);
        } else {
            command = String::from(arg);  // merge, compare, rebuild, compare
        }
    }

    println!("Going to load {} trees.", trees.len());

    if command == String::from("merge") {
        let distance_matrices: Vec<TreeDistanceMatrix>;

        if mapping.len() == 0 {
            println!("No mapping specified.");
            distance_matrices = trees.par_iter().map(
                |t| global_supertrees::utils::convert_file_to_distance_matrix(t.to_string())
            ).collect();
        } else {
            distance_matrices = trees.par_iter().map(
                |t| global_supertrees::utils::convert_file_to_distance_matrix(t.to_string()).merge_organisms(&mapping)
            ).collect();
        }


//        for dm in distance_matrices.iter() {
//            let mean = dm.distance_matrix.sum() / (dm.shape()[0] as f32 * dm.shape()[0] as f32);
//            let mut min = f32::MAX;
//            let mut max = 0.0;
//
//            for el in dm.distance_matrix.iter() {
//                if el < &min {
//                    min = *el;
//                }
//                if el > &max {
//                    max = *el;
//                }
//            }
//
//            let mut values : Vec<f32> = dm.distance_matrix.iter().map(|x| *x).collect();
//            values.sort_by(|a, b| a.partial_cmp(b).unwrap());
//            let median = values[values.len() / 2];
//            println!("Mean: {}, max: {}, min: {}, Median: {}", mean, max, min, median);
//        }
//
//        std::process::exit(0);

        let merged_tree = mean_merge_distance_matrices(distance_matrices);

        let mut out_tree_file = File::create("merged.tree").expect("IO Error while creating file.");
        out_tree_file.write_all(merged_tree.to_newick().as_bytes()).expect("IO Error while writing merged tree.");
        println!("Created tree with {} leaves.", merged_tree.get_leaves().len());

        let final_distance_matrix = merged_tree.to_distance_matrix();
        let mut out_dm_file = File::create("merged_distance_matrix.csv").expect("IO Error while creating file.");

        out_dm_file.write_all(final_distance_matrix.to_csv().as_bytes()).expect("IO Error while writing merged distance matrix.");
    } else if command == String::from("rebuild") {
        let rebuilt_trees: Vec<Tree> = trees.par_iter().map(
            |t| global_supertrees::utils::convert_file_to_distance_matrix(t.to_string()).to_tree()
        ).collect();
        for (t, tree) in rebuilt_trees.iter().enumerate() {
            let mut out_tree_file = File::create(format!("rebuilt_tree_{}.tree", t)).expect("IO Error while creating file.");
            out_tree_file.write_all(tree.to_newick().as_bytes());
        }
    } else if command == String::from("compare") {
        println!("Comparing");
        if trees.len() == 2 {
            println!("Converting");
            let distance_matrices : Vec<TreeDistanceMatrix> = trees.par_iter().map(
                |t| global_supertrees::utils::convert_file_to_distance_matrix(t.to_string())
            ).collect();


            let mut out_dm_file_a = File::create("dm_a.csv").expect("IO Error while creating file.");
            out_dm_file_a.write_all(distance_matrices[0].to_csv().as_bytes()).expect("IO Error while writing merged distance matrix.");

            let mut out_dm_file_b = File::create("dm_b.csv").expect("IO Error while creating file.");
            out_dm_file_b.write_all(distance_matrices[1].to_csv().as_bytes()).expect("IO Error while writing merged distance matrix.");

        } else {
            println!("Too many trees, 2 are needed.");
        }
    } else {
        println!("Invalid command: {}", command);
    }

}
