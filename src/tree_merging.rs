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

    sibling_merging(new_trees);

}

fn get_tree_size(tree : &Tree) -> usize {
    return tree.traverse_children().len();
}


fn sibling_merging(mut trees : Vec<Tree>) -> () {
    let max_ancestor_search = 1;

    trees.sort_unstable_by(
        |a , b| get_tree_size(a).cmp(&get_tree_size(b))
    );
    let mut base_tree : Tree;

    match trees.pop() {
        Some(t) => { base_tree = t; },
        None => { return; } // TODO: return an empty tree;
    }
    println!("{:#?}", base_tree);

    trees.reverse();  // Start with the largest trees.

    for _anc in 0..max_ancestor_search {
        for tree in trees.iter() {
            //println!("{:?}", tree);
            for (child, _dist) in tree.traverse_children() {
                let leaf_set = utils::vec_to_set(&child.leaves);
                for (bi, (base_child, _base_dist)) in base_tree.traverse_children_snapshot().iter().enumerate() {
                    let base_leaf_set = utils::vec_to_set(&base_child.leaves);
                    //println!("base_leaf_set: {:?}", base_leaf_set);
                    let overlapping_leaves : Vec<&String> = leaf_set.intersection(&base_leaf_set).collect();
                    if overlapping_leaves.len() == 0 {
                        continue;  // Skip if there are no overlapping siblings.
                    }
                    let new_leaves : Vec<&String> = leaf_set.difference(&base_leaf_set).collect();
                    if new_leaves.len() == 0 {
                        continue; // Skip if there are no new leaves.
                    }

                    //println!("Found new siblings! {:?} in {:?}, {:?}", new_leaves, base_leaf_set, leaf_set);

                    for new_leaf in new_leaves.iter() {
                        match child.leaves.iter().position(| c | &c == new_leaf) {
                            Some(idx) => {
                                let new_leaf_distance = child.leaf_distances[idx];
                                base_tree.perform_operation_on_branch(
                                    bi,
                                    0,
                                    &mut | t : &mut Tree | t.add_leaf(new_leaf.to_string(), new_leaf_distance)
                                );
                            }
                            None => { continue; }
                        }
                    }

                }
            }
        }
    }
    println!("{:#?}", base_tree);
}


fn merge_distance_matrices(dms : Vec<TreeDistanceMatrix>) -> TreeDistanceMatrix {
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