use rayon::prelude::*;
use ndarray::prelude::*;
use std::collections::HashSet;
use std::f32;
use crate::graph_tree::*;
use crate::tree_distance_matrix::*;
use crate::utils;
use rand::prelude::*;
use rand::seq::SliceRandom;
use petgraph::visit::Dfs;
use petgraph::graph::node_index;

pub fn merge_trees(trees: Vec<Tree>) -> Tree {
    let distance_matrices: Vec<TreeDistanceMatrix> =
        trees.par_iter().map(|t| t.to_distance_matrix()).collect();

    mean_merge_distance_matrices(distance_matrices)
}

pub fn mean_merge_distance_matrices(distance_matrices : Vec<TreeDistanceMatrix>) -> Tree {
    println!("Started mean merging");
    let leaves: Vec<Vec<String>> = distance_matrices
        .par_iter()
        .map(|t| utils::keys_vec(&t.leaf_map))
        .collect();

    println!("Merging distance matrices");
    // Create a shared matrix from which every tree can get the new mean distances.
    let merged_distances = merge_distance_matrices(distance_matrices);

    println!("Replacing values with mean values.");
    let new_trees: Vec<Tree> = leaves
        .par_iter()
        .map(|t| {
            merged_distances
                .get_partial_distance_matrix(t)
                .to_tree()
        })
        .collect();

//    let new_trees : Vec<Tree> = distance_matrices
//        .par_iter()
//        .map(|t| {
//            t.neighbour_joining()
//        }).collect();

    println!("Started sibling merging");
    let mut merged_trees : Vec<Tree> = vec![0, 1, 2, 3, 4].par_iter().map(
        |seed| sibling_merging(new_trees.clone(), 3, *seed)
    ).collect();
    println!();

    merged_trees.sort_by(
        |a, b| a.get_leaves().len().cmp(&b.get_leaves().len())
    );

    for tree in merged_trees.iter() {
        print!("{}, ", tree.get_leaves().len());
    }

    merged_trees.pop().unwrap()
}


fn sibling_merging(mut trees: Vec<Tree>, bootstrap_n: usize, seed: usize) -> Tree {
    let max_ancestor_search = 3;

    let mut base_tree : Tree;
    let mut tree : &Tree;
    let n_trees = trees.len();
    let n_iter = n_trees * bootstrap_n;

    let seed_arr: [u8; 32] = [seed as u8; 32];
    let mut rng: StdRng = SeedableRng::from_seed(seed_arr);

    trees.shuffle(&mut rng);

    match trees.pop() {
        Some(t) => { base_tree = t; }
        None => { return Tree::new(); }  // Return an empty tree if the list is empty
    }

    let mut unmatched_leaves : HashSet<String> = HashSet::new();
    for tree in trees.iter() {
        for leaf in tree.get_leaves() {
            unmatched_leaves.insert(leaf);
        }
    }

    base_tree.to_distance_matrix().distance_matrix.map(|x : &f32| assert!(!x.is_nan()));

    for anc in 1..(max_ancestor_search + 1) {  // Iterates up the levels in the tree where siblings are shared
        let mut iteration_added_siblings = true;
        while iteration_added_siblings {  // Iteratively attempt to add more siblings
            iteration_added_siblings = false;
            println!("==({})==", anc);

            for _ in 0..n_iter {
                tree = trees.iter().choose(&mut rng).unwrap();
                let mut dfs = Dfs::new(&tree.graph,  node_index(0));

                while let Some(node) = dfs.next(&tree.graph) {
                    let base_leaves = base_tree.get_leaves();

                    let node_name = tree.graph.node_weight(node).unwrap();
                    if node_name == "<root>" {
                        continue;
                    }

                    if node_name.starts_with(">>") {
                        // This is a parent node.
                    } else {
                        // Encountered a leaf node.
                        if !base_leaves.contains(node_name) {
                            //println!("No overlap between base tree and {}", node_name);
                            continue;  // Skip if there is no overlap.
                        }
                        unmatched_leaves.remove(node_name);

                        let base_node = base_tree.find_node_idx(node_name).unwrap();
                        let siblings = Tree::get_node_siblings(&tree.graph, node, anc);

                        for sibling in siblings {
                            let sibling_name = tree.graph.node_weight(sibling).unwrap();
                            if !base_leaves.contains(sibling_name) {
                                let sibling_distance = Tree::get_node_distance_to_parent(&tree.graph,sibling, anc);
                                base_tree.add_sibling_n_removed(
                                    base_node,
                                    sibling_name.to_string(),
                                    sibling_distance,
                                    anc
                                );
                                //println!("Added sibling {}", sibling_name);
                                unmatched_leaves.remove(sibling_name);
                                iteration_added_siblings = true;
                            }
                        }
                    }
                }
            }
        }
    }

    base_tree.to_distance_matrix().distance_matrix.map(|x : &f32| assert!(!x.is_nan()));

    println!("Unmatched leaves: {:?}", unmatched_leaves);

    //println!("{:#?}", base_tree);
    return base_tree;
}

#[allow(dead_code)]
/// This function creates a merged distance matrix, where some fields are filled with NaN because
/// there is no mapping between the species in any of the trees. Trees can be reconstructed using
/// the `get_partial_distance_matrix` function.
fn merge_distance_matrices(dms: Vec<TreeDistanceMatrix>) -> TreeDistanceMatrix {
    let mut all_leaves: HashSet<String> = HashSet::new();
    for distance_matrix in dms.iter() {
        for leaf in distance_matrix.leaf_map.keys() {
            all_leaves.insert(leaf.clone());
        }
    }

    let mut overlapping_distance_matrix: Array2<f32> =
        Array2::zeros((all_leaves.len(), all_leaves.len()));
    overlapping_distance_matrix.fill(f32::NAN);
    let mut count_matrix: Array2<f32> = Array2::zeros((all_leaves.len(), all_leaves.len()));

    for distance_matrix in dms.iter() {
        for (x, leaf_x) in all_leaves.iter().enumerate() {
            for (y, leaf_y) in all_leaves.iter().enumerate() {
                let dist = distance_matrix.get_distance_ref(leaf_x, leaf_y);
                match dist {
                    Some(d) => {
                        overlapping_distance_matrix[[x, y]] = d;
                        count_matrix[[x, y]] += 1.0
                    }
                    None => {}
                }
            }
        }
    }

    // Get the average distances.
    let average_matrix = overlapping_distance_matrix / count_matrix;

    return TreeDistanceMatrix::new_from_matrix_and_leaves(
        average_matrix,
        all_leaves.iter().map(|l| l.clone()).collect(),
    );
}
