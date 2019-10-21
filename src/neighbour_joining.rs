use ndarray::prelude::*;
use std::f32;
use petgraph::graph::node_index;
use petgraph::Graph;
use petgraph::prelude::NodeIndex;
use crate::tree_distance_matrix::TreeDistanceMatrix;
use crate::graph_tree::Tree;


/// Calculates a Q matrix based on a distance matrix for the neighbour joining algorithm.
fn calculate_q_matrix(dm: &Array2<f32>) -> Array2<f32> {
    let mut q_matrix: Array2<f32> = Array2::zeros(dm.raw_dim());
    let n = dm.shape()[0]; // Number of taxa

    let row_sums : Vec<f32> = dm.sum_axis(Axis(1)).to_vec();
    let col_sums : Vec<f32> = dm.sum_axis(Axis(0)).to_vec();
    let mult_dm = dm * (n as f32 - 2.0);

    for i in 0..n {
        for j in 0..n {
            q_matrix[[i, j]] = mult_dm[[i, j]] - row_sums[i] - col_sums[j];
        }
    }
    q_matrix
}

/// Calculates a pair distance matrix based on a Qmatrix for the neighbour joining algorithm.
/// Results in a new distance matrix with the new pair as a single node.
fn calculate_pair_distance_matrix(i: usize, j: usize, dm: &Array2<f32>) -> Array2<f32> {
    assert_eq!(dm.shape()[0] > 0, true);
    let mut new_distances: Array2<f32> = Array2::zeros((dm.shape()[0] - 1, dm.shape()[1] - 1));
    // New indexes
    let mut x_n = 0;
    let mut y_n;
    let final_idx = new_distances.shape()[0] - 1;

    for x in 0..dm.shape()[0] {
        y_n = 0;
        for y in 0..dm.shape()[1] {
            if (x == i) | (x == j) {
                new_distances[[final_idx, y_n]] = 0.5 * (dm[[i, y]] + dm[[j, y]] - dm[[i, j]]);
                if (y != i) & (y != j) {
                    y_n += 1;
                }
            } else if (y == i) | (y == j) {
                new_distances[[x_n, final_idx]] = 0.5 * (dm[[i, x]] + dm[[j, x]] - dm[[i, j]]);
            } else {
                new_distances[[x_n, y_n]] = dm[[x, y]];
                y_n += 1;
            }
        }
        if (x != i) & (x != j) {
            x_n += 1;
        }
    }

    return new_distances;
}

/// Finds the leaf distances for the pair inside the Tree node.
fn calculate_new_leaf_distances(f: usize, g: usize, dm: &Array2<f32>) -> (f32, f32) {
    // f and g are the paired taxa and u is the new created node.
    let n = dm.shape()[0] as f32;
    let sum_fk = dm.row(f).sum();
    let sum_gk = dm.row(g).sum();
    let dist_fg = dm[[f, g]];

    let sum_block = sum_fk - sum_gk;
    let base_block = 0.5 * dist_fg;
    let over_block = 1.0 / (2.0 * (n - 2.0));

    let mut d_fu = base_block + (over_block * sum_block);
    if d_fu.is_nan() {
        return (0.0, 0.0);
    }

    let mut d_gu = dist_fg - d_fu;

    // As per http://www.icp.ucl.ac.be/~opperd/private/neighbor.html
    // Negative distances are possible in neighbour joining. You can correct for this by setting
    // the distance to zero and adding it to the adjacent branch.
    if d_fu < 0.0 {
        d_gu += d_fu * -1.0;
        d_fu = 0.0;
    } else if d_gu < 0.0 {
        d_fu += d_gu * -1.0;
        d_gu = 0.0;
    }

    return (d_fu, d_gu);
}

/// Returns the lowest index in the matrix where i != j
fn find_min_matrix_ne(matrix: Array2<f32>) -> [usize; 2] {

    let mut lowest: f32 = f32::MAX;
    let mut idx = [0, 0];
    let mut v: f32;
    for i in 0..matrix.shape()[0] {
        for j in 0..matrix.shape()[1] {
            if i != j {
                v = matrix[[i, j]];
                if v < lowest {
                    lowest = v;
                    idx = [i, j];
                }
            }
        }
    }
    return idx;
}

/// Performs the neighbour joining algorithm on this distance matrix to create a tree.
pub fn neighbour_joining(dm : &TreeDistanceMatrix) -> Tree {
    let mut distance_matrix = dm.distance_matrix.clone();
    println!("Started neighbour joining");
    let mut leaf_map_inv_vec: Vec<(usize, String)> = dm
        .leaf_map_inv
        .iter()
        .map(|val| (*val.0, val.1.clone()))
        .collect();

    leaf_map_inv_vec.sort_by(|a, b| a.0.cmp(&b.0));

    // This will contain all the newly created branches
    let mut branches: Vec<NodeIndex> = Vec::new();

    let mut graph : Graph<String, f32> = Graph::new();
    // Add root node on index 0
    graph.add_node(String::from("<root>"));

    for (_, leaf) in leaf_map_inv_vec.iter() {
        // Add nodes in the proper order, their NodeIndex will be their leaf index + 1
        graph.add_node(leaf.to_string());
    }
    let mut n_iters = 1;

    while distance_matrix.shape()[0] > 1 {
        //println!("{:?} - {}", distance_matrix.shape(), branches.len());

        // Calculate the Q matrix
        let q_matrix = calculate_q_matrix(&distance_matrix);
        // Find the pair of leaves that have the lowest value in the q_matrix
        let lowest_pair = find_min_matrix_ne(q_matrix);

        // Calculate the distance matrix from each of the taxa to the new node
        // In this pair_distance_matrix, the final row/col is now designated the previously
        // selected pair.
        let pair_distance_matrix = calculate_pair_distance_matrix(
            lowest_pair[0],
            lowest_pair[1],
            &distance_matrix,
        );

        // Calculate the distances between the new node and the two merged taxa.
        let pair_leaf_distances = calculate_new_leaf_distances(
            lowest_pair[0],
            lowest_pair[1],
            &distance_matrix,
        );

        let mut children : Vec<NodeIndex> = Vec::new();
        let mut child_distances : Vec<f32> = Vec::new();
        let mut rem_leaves: Vec<usize> = Vec::new();
        let mut rem_branches: Vec<usize> = Vec::new();

        // Check if the first node is already a branch
        if lowest_pair[0] >= leaf_map_inv_vec.len() {
            let branch_index = lowest_pair[0] - leaf_map_inv_vec.len();
            rem_branches.push(branch_index);
            children.push(branches[branch_index]);
            child_distances.push(pair_leaf_distances.0);
        } else {
            rem_leaves.push(lowest_pair[0]);
            children.push(node_index(leaf_map_inv_vec.get(lowest_pair[0]).expect(
                "Node with first lowest pair ID does not exist.").0 + 1));
            child_distances.push(pair_leaf_distances.0);
        }

        // Check if the second node is already a tree
        if lowest_pair[1] >= leaf_map_inv_vec.len() {
            let branch_index = lowest_pair[1] - leaf_map_inv_vec.len();
            rem_branches.push(branch_index);
            children.push(branches[branch_index]);
            child_distances.push(pair_leaf_distances.1);
        } else {
            rem_leaves.push(lowest_pair[1]);
            children.push(node_index(leaf_map_inv_vec.get(lowest_pair[1]).expect(
                "Node with second lowest pair ID does not exist.").0 + 1));
            child_distances.push(pair_leaf_distances.1);
        }

        // Remove leaves from the leaf map if they have been selected
        rem_leaves.sort_unstable_by(|a, b| b.cmp(a));
        for ri in rem_leaves {
            leaf_map_inv_vec.remove(ri);
        }

        // Remove branches that have been selected.
        rem_branches.sort_unstable_by(|a, b| b.cmp(a));
        for ri in rem_branches {
            branches.remove(ri);
        }

        assert_eq!(children.len(), 2);
        assert_eq!(child_distances.len(), 2);

        // Create connections
        let new_branch = graph.add_node(format!(">>{}", n_iters));
        graph.add_edge(new_branch, children[0], child_distances[0]);
        graph.add_edge(new_branch, children[1], child_distances[1]);

        distance_matrix = pair_distance_matrix;

        branches.push(new_branch);

        n_iters += 1;
    }

    println!("Finished neighbour joining.");
    if branches.len() >= 1 {
        graph.add_edge(node_index(0), branches[(branches.len() - 1)], f32::NAN);
    }

    Tree { graph }
}
