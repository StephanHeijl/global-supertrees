use std::f32;
use std::collections::HashMap;
use std::collections::HashSet;
use ndarray::prelude::*;
use tree::*;
use tree_distance_matrix::*;
use rayon::prelude::*;
use::utils;


pub fn merge_trees(trees : Vec<Tree>) {
    let distance_matrices : Vec<TreeDistanceMatrix> = trees.iter().map(| t | t.to_distance_matrix()).collect();

    let tree_leafs : Vec<Vec<String>> = distance_matrices.iter().map(
        |t| utils::keys_vec(&t.leaf_map)
    ).collect();
    // Create a shared matrix from which every tree can get the new mean distances.
    let merged_distances = merge_distance_matrices(distance_matrices);

    let new_trees : Vec<Tree> = tree_leafs.par_iter().map(
        | t | merged_distances.get_partial_distance_matrix(t).neighbour_joining()
    ).collect();

}

pub fn merge_distance_matrices(dms : Vec<TreeDistanceMatrix>) -> TreeDistanceMatrix {
    let mut all_leaves : HashSet<String> = HashSet::new();
    for distance_matrix in dms.iter() {
        for leaf in distance_matrix.leaf_map.keys() {
            all_leaves.insert(leaf.clone());
        }
    }

    let mut overlapping_distance_matrix : Array2<f32> = Array2::zeros((all_leaves.len(), all_leaves.len()));
    overlapping_distance_matrix.fill(f32::NAN);
    let mut count_matrix : Array2<f32> = Array2::zeros((all_leaves.len(), all_leaves.len()));

    for distance_matrix in dms.iter() {
        for (x, leaf_x) in all_leaves.iter().enumerate() {
            for (y, leaf_y) in all_leaves.iter().enumerate() {
                let dist = distance_matrix.get_distance_ref(leaf_x, leaf_y);
                match dist {
                    Some(d) => { overlapping_distance_matrix[[x, y]] = d; count_matrix[[x, y]] += 1.0 },
                    None => { }
                }
            }
        }
    }

    // Get the average distances.
    let average_matrix = overlapping_distance_matrix / count_matrix;
    return TreeDistanceMatrix::new_from_matrix_and_leaves(
        average_matrix,
        all_leaves.iter().map(|l| l.clone()).collect()
    );

}