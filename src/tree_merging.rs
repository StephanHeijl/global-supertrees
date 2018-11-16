use std::f32;
use std::collections::HashMap;
use std::collections::HashSet;
use ndarray::prelude::*;
use tree::*;
use tree_distance_matrix::*;
use rayon::prelude::*;


pub fn merge_trees(trees : Vec<Tree>) {
    let distance_matrices : Vec<TreeDistanceMatrix> = trees.iter().map(| t | t.to_distance_matrix()).collect();

    merge_distance_matrices(distance_matrices);
}

pub fn merge_distance_matrices(dms : Vec<TreeDistanceMatrix>) {
    let mut all_leaves : HashSet<String> = HashSet::new();
    for distance_matrix in dms.iter() {
        for leaf in distance_matrix.leaf_map.keys() {
            all_leaves.insert(leaf.clone());
        }
    }

    println!("{:?}", all_leaves);

    let mut overlapping_distance_matrix : Array2<f32> = Array2::zeros((all_leaves.len(), all_leaves.len()));

    println!("{:?}", overlapping_distance_matrix.shape());

    for distance_matrix in dms.iter() {
        for (x, leaf_x) in all_leaves.iter().enumerate() {
            for (y, leaf_y) in all_leaves.iter().enumerate() {
                distance_matrix.get_distance_ref(leaf_x, leaf_y);
            }
        }
    }

}