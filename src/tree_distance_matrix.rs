//! Tree Distance Matrix Module
use ndarray::prelude::*;
use ndarray::stack;
use rayon::prelude::*;
use std::cmp::min;
use std::collections::HashMap;
use std::collections::BTreeMap;
use std::f32;

use regex::Regex;

use petgraph::graph::node_index;
use petgraph::Graph;
use petgraph::prelude::NodeIndex;

use graph_tree::*;
use uniprot::*;
use utils;

#[derive(Debug)]
pub struct TreeDistanceMatrix {
    pub leaf_map: HashMap<String, usize>,
    pub leaf_map_inv: HashMap<usize, String>,
    pub distance_matrix: Array2<f32>,
}

/// Tree distance matrix processing implementation
impl TreeDistanceMatrix {

    /// Calculates a Q matrix based on a distance matrix for the neighbour joining algorithm.
    pub fn calculate_q_matrix(dm: &Array2<f32>) -> Array2<f32> {
        let mut q_matrix: Array2<f32> = Array2::zeros(dm.raw_dim());
        let n = dm.shape()[0]; // Number of taxa

        for i in 0..n {
            for j in 0..n {
                q_matrix[[i, j]] = ((n as f32 - 2.0) * dm[[i, j]])
                    - dm.row(i).sum()
                    - dm.column(j).sum();
            }
        }
        q_matrix
    }

    /// Calculates a pair distance matrix based on a Qmatrix for the neighbour joining algorithm.
    /// Results in a new distance matrix with the new pair as a single node.
    pub fn calculate_pair_distance_matrix(i: usize, j: usize, dm: &Array2<f32>) -> Array2<f32> {
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
    pub fn calculate_new_leaf_distances(f: usize, g: usize, dm: &Array2<f32>) -> (f32, f32) {
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
    pub fn find_min_matrix_ne(matrix: Array2<f32>) -> [usize; 2] {

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
    pub fn neighbour_joining(&self) -> Tree {
        println!("Start neighbour joining");
        let mut distance_matrix = self.distance_matrix.clone();
        let mut leaf_map_inv_vec: Vec<(usize, String)> = self
            .leaf_map_inv
            .iter()
            .map(|val| (*val.0, val.1.clone()))
            .collect();

        leaf_map_inv_vec.sort_by(|a, b| a.0.cmp(&b.0));

        // This will contain all the newly created branches
        let mut branches: HashMap<usize, NodeIndex> = HashMap::new();

        let mut graph : Graph<String, f32> = Graph::new();
        // Add root node on index 0
        graph.add_node(String::from("<root>"));

        for (_, leaf) in leaf_map_inv_vec.iter() {
            // Add nodes in the proper order, their NodeIndex will be their leaf index + 1
            graph.add_node(leaf.to_string());
        }
        let mut n_iters = 1;

        while distance_matrix.shape()[0] > 1 {
            println!("{}", n_iters);
            // Calculate the Q matrix
            let q_matrix = TreeDistanceMatrix::calculate_q_matrix(&distance_matrix);
            // Find the pair of leaves that have the lowest value in the q_matrix
            let lowest_pair = TreeDistanceMatrix::find_min_matrix_ne(q_matrix);

            // Calculate the distance matrix from each of the taxa to the new node
            // In this pair_distance_matrix, the final row/col is now designated the previously
            // selected pair.
            let pair_distance_matrix = TreeDistanceMatrix::calculate_pair_distance_matrix(
                lowest_pair[0],
                lowest_pair[1],
                &distance_matrix,
            );

            // Calculate the distances between the new node and the two merged taxa.
            let pair_leaf_distances = TreeDistanceMatrix::calculate_new_leaf_distances(
                lowest_pair[0],
                lowest_pair[1],
                &distance_matrix,
            );

            let mut children : Vec<NodeIndex> = Vec::new();
            let mut child_distances : Vec<f32> = Vec::new();
            let mut rem_leaves: Vec<usize> = Vec::new();

            // Check if the first node is already a branch
            if branches.contains_key(&lowest_pair[0]) {
                children.push(branches.remove(&lowest_pair[0]).expect("Branch with first lowest pair ID does not exist."));
                child_distances.push(pair_leaf_distances.0);
            } else {
                rem_leaves.push(lowest_pair[0]); 
                children.push(node_index(leaf_map_inv_vec.get(lowest_pair[0]).expect("Node with first lowest pair ID does not exist.").0 + 1));
                child_distances.push(pair_leaf_distances.0);
            }

            // Check if the second node is already a tree
            if branches.contains_key(&lowest_pair[1]) {
                children.push(branches.remove(&lowest_pair[1]).expect("Branch with second lowest pair ID does not exist."));
                child_distances.push(pair_leaf_distances.1);
            } else {
                rem_leaves.push(lowest_pair[1]);
                children.push(node_index(leaf_map_inv_vec.get(lowest_pair[1]).expect("Node with second lowest pair ID does not exist.").0 + 1));
                child_distances.push(pair_leaf_distances.1);
            }

            // Remove leaves from the leaf map if they have been selected
            let mut leaf_map_inv_vec_new = Vec::new();
            for (i, el) in leaf_map_inv_vec.clone().iter().enumerate() {
                if !rem_leaves.contains(&i) {
                    leaf_map_inv_vec_new.push(el.clone());
                }
            }
            leaf_map_inv_vec = leaf_map_inv_vec_new.clone();

            // Create connections
            let new_branch = graph.add_node(format!(">>{}", n_iters));
            graph.add_edge(new_branch, children[0], child_distances[0]);
            graph.add_edge(new_branch, children[1], child_distances[1]);

            // Shift all the branches to conform to the new shifted index
            // as the distance matrix has shrunk
            let mut branches_new: HashMap<usize, NodeIndex> = HashMap::new();
            // Get the tree keys without borrowing them.
            let branch_keys: Vec<usize> = branches.keys().map(|x| *x).collect();
            for idx in branch_keys {
                branches_new.insert(idx - 2, branches.remove(&idx).expect("Could not move branch"));
            }
            branches = branches_new;

            distance_matrix = pair_distance_matrix;

            branches.insert(distance_matrix.shape()[0] - 1, new_branch);

            n_iters += 1;
        }

        graph.add_edge(node_index(0), branches[&(branches.len() - 1)], f32::NAN);

        Tree { graph }
    }

    /// Returns the distance between two leaves.
    #[allow(dead_code)]
    pub fn get_distance(&self, leaf_one: String, leaf_two: String) -> Option<f32> {
        return self.get_distance_ref(&leaf_one, &leaf_two);
    }

    /// Returns the distance between two leaves based on references.
    #[allow(dead_code)]
    pub fn get_distance_ref(&self, leaf_one: &String, leaf_two: &String) -> Option<f32> {
        let l1: Option<&usize> = self.leaf_map.get(leaf_one);
        let l2: Option<&usize> = self.leaf_map.get(leaf_two);
        let x: usize;
        let y: usize;
        match l1 {
            Some(i) => x = *i,
            None => {
                return None;
            }
        }
        match l2 {
            Some(i) => y = *i,
            None => {
                return None;
            }
        }
        return Some(self.distance_matrix[[x, y]]);
    }

    /// Converts the distance matrix to a CSV string.
    #[allow(dead_code)]
    pub fn to_csv(&self) -> String {
        let mut csv = String::from("");
        let delim = ',';

        let mut max_key: usize = 0;
        for key in self.leaf_map_inv.keys() {
            if key > &max_key {
                max_key = *key;
            }
        }

        let mut leaves: Vec<String> = Vec::new();

        csv.push_str(",");
        for l in 0..max_key + 1 {
            match self.leaf_map_inv.get(&l) {
                Some(name) => {
                    csv.push_str(name);
                    leaves.push(name.to_string());
                    if l < max_key {
                        csv.push(delim);
                    }
                }
                None => {}
            }
        }
        csv.push('\n');

        for (r, row) in self.distance_matrix.outer_iter().enumerate() {
            csv.push_str(&leaves[r]);
            csv.push(delim);
            let mut v = 0;
            for val in row {
                csv.push_str(&val.to_string());
                if v < leaves.len() - 1 {
                    csv.push(delim);
                }
                v += 1;
            }
            csv.push('\n');
        }
        return csv;
    }

    /// Gets the shape of the distance matrix.
    #[allow(dead_code)]
    pub fn shape(&self) -> &[usize] {
        return self.distance_matrix.shape();
    }


    /// Returns the index of the first zero in a 1D ArrayView from position 1.
    /// Returns 0 if no zero has been found.
    #[allow(dead_code)]
    pub fn find_first_zero(identity_row: ArrayView1<usize>) -> usize {
        let mut i = 0;
        for n in identity_row.iter() {
            if i == 0 {
                i += 1;
                continue;
            }
            if *n == 0 {
                return i;
            }
            i += 1;
        }
        return 0;
    }

    /// Returns the indices of the last zeros for two rows in the identity matrix based on the 2 leaf identifiers (l1, l2).
    fn find_final_parent(l1: usize, l2: usize, identity_matrix: &Array2<usize>) -> (usize, usize) {
        let p1 = TreeDistanceMatrix::find_first_zero(identity_matrix.slice(s![l1, ..]));
        let p2 = TreeDistanceMatrix::find_first_zero(identity_matrix.slice(s![l2, ..]));

        (p1, p2)
    }

    ///  Find the first common ancestor of two leaves.  If they do not share any explicit ancestors,
    ///  the first common ancestor becomes the root of the tree (0). If the leaves are on the same (sub)tree,
    ///  their shared subtree is returned. The ancestor trees are returned as a number representing the subtree id.
    pub fn find_first_common_ancestor(
        l1: usize,
        l2: usize,
        identity_matrix: &Array2<usize>,
    ) -> usize {

        let id_row_1 = identity_matrix.slice(s![l1, ..]);
        let id_row_2 = identity_matrix.slice(s![l2, ..]);

        if id_row_1 == id_row_2 {
            return TreeDistanceMatrix::find_final_parent(l1, l2, identity_matrix).0 - 1;
        }

        for i in 1..id_row_1.len() {
            if id_row_1[i] != id_row_2[i] {
                return i - 1;
            }
        }
        return 0;
    }

    /// Generates a part of a distance matrix for use with multicore processing.
    #[cfg(not(feature = "singlecore"))]
    fn generate_partial_distance_matrix(
        leaf_distance_matrix: &Array2<f32>,
        identity_matrix: &Array2<usize>,
        iteration: usize,
        size: usize,
    ) -> Array2<f32> {

        let mut distance_matrix: Array2<f32> = Array2::zeros((size, size));
        let n_leaves = identity_matrix.shape()[0];

        let n_partials_axis = ((n_leaves as f32) / (size as f32)).ceil() as usize;

        let x_iter = iteration % n_partials_axis;
        let y_iter = iteration / n_partials_axis;

        let x_start = x_iter * size;
        let x_stop = min(x_iter * size + size, n_leaves);
        let y_start = y_iter * size;
        let y_stop = min(y_iter * size + size, n_leaves);

        for x in x_start..x_stop {
            for y in y_start..y_stop {
                let xd = x - x_start;
                let yd = y - y_start;

                if x == y {
                    distance_matrix[[xd, yd]] = 0.0;
                    continue;
                }

                let (xi, yi) = TreeDistanceMatrix::find_final_parent(x, y, identity_matrix);
                let fca = TreeDistanceMatrix::find_first_common_ancestor(x, y, identity_matrix);

                // Find the total distance from each leaf to the root.
                let leaf_root_distance: f32 =
                    leaf_distance_matrix[[x, xi]] + leaf_distance_matrix[[y, yi]];

                // Find the distance between the root and the first common ancestor.
                let root_fca_distance: f32 = leaf_distance_matrix[[x, fca]];
                let distance = f32::abs(leaf_root_distance - (root_fca_distance * 2.0));

                // Add the distance to the distance matrix
                distance_matrix[[xd, yd]] = distance;
            }
        }

        distance_matrix
    }

    /// Generates a full distance matrix in multiprocessing mode with rayon.
    #[cfg(not(feature = "singlecore"))]
    fn generate_full_distance_matrix(
        leaf_distance_matrix: Array2<f32>,
        identity_matrix: Array2<usize>,
        n_leaves: usize,
    ) -> Array2<f32> {
        /* Generates a full distance matrix */

        //let mut distance_matrices : Vec<Array2<f32>> = Vec::new();
        let max_size = 1024 as f32;
        let parts = ((n_leaves as f32) / max_size).ceil() as usize;

        let starts: Vec<usize> = (0..parts.pow(2)).collect();

        let distance_matrices: Vec<Array2<f32>> = starts
            .par_iter()
            .map(|&s| {
                TreeDistanceMatrix::generate_partial_distance_matrix(
                    &leaf_distance_matrix,
                    &identity_matrix,
                    s,
                    max_size as usize,
                )
            })
            .collect();

        // Compile partial distance matrices
        let mut rows: Vec<Array2<f32>> = Vec::new();
        for x in 0..parts {
            let start_x = x * parts;
            let stop_x = x * parts + parts;
            let row: Vec<ndarray::ArrayBase<ndarray::ViewRepr<&f32>, ndarray::Dim<[usize; 2]>>> =
                distance_matrices[start_x..stop_x]
                    .iter()
                    .map(|dm| dm.view())
                    .collect();
            rows.push(stack(Axis(1), &row).expect("Stacking on axis 1 failed."));
        }

        let row_views: Vec<ndarray::ArrayBase<ndarray::ViewRepr<&f32>, ndarray::Dim<[usize; 2]>>> =
            rows.iter().map(|row| row.view()).collect();

        let distance_matrix = stack(Axis(0), &row_views).expect("Stacking on axis 0 failed.");
        distance_matrix.slice(s![..n_leaves, ..n_leaves]).to_owned()
    }

    /// Generates a full distance matrix in singlecore mode.
    #[cfg(feature = "singlecore")]
    fn generate_full_distance_matrix(
        leaf_distance_matrix: Array2<f32>,
        identity_matrix: Array2<usize>,
        n_leaves: usize,
    ) -> Array2<f32> {
        /* Generates a full distance matrix */

        let mut distance_matrix: Array2<f32> = Array2::zeros((n_leaves, n_leaves));

        for x in 0..n_leaves {
            /* Iterate over each leaf for the X axis. */
            for y in 0..n_leaves {
                /* Iterate over each leaf for the Y axis */
                if x == y {
                    distance_matrix[[x, y]] = 0.0;
                    continue;
                }

                let (xi, yi) = TreeDistanceMatrix::find_final_parent(x, y, &identity_matrix);
                let fca = TreeDistanceMatrix::find_first_common_ancestor(x, y, &identity_matrix);

                // Find the total distance from each leaf to the root.
                let leaf_root_distance: f32 =
                    leaf_distance_matrix[[x, xi]] + leaf_distance_matrix[[y, yi]];

                // Find the distance between the root and the first common ancestor.
                let root_fca_distance: f32 = leaf_distance_matrix[[x, fca]];
                let distance = f32::abs(leaf_root_distance - (root_fca_distance * 2.0));
                distance_matrix[[x, y]] = distance;
            }
        }

        distance_matrix
    }

    /// Returns a part of the distance matrix including only the leaves specified.
    pub fn get_partial_distance_matrix(&self, leaves: &Vec<String>) -> TreeDistanceMatrix {
        let mut leaf_map: HashMap<String, usize> = HashMap::new();
        let mut leaf_map_inv: HashMap<usize, String> = HashMap::new();
        let self_leaves: Vec<String> = utils::keys_vec(&self.leaf_map);

        let mut l = 0;
        let mut dm_idx: Vec<usize> = Vec::new();
        for leaf in leaves.iter() {
            if self_leaves.contains(&leaf) {
                leaf_map.insert(leaf.clone(), l);
                leaf_map_inv.insert(l, leaf.clone());
                dm_idx.push(self.leaf_map[leaf]);
                l += 1;
            }
        }

        let mut ndm: Array2<f32> = Array2::zeros((l, l));

        for (x, ix) in dm_idx.iter().enumerate() {
            for (y, iy) in dm_idx.iter().enumerate() {
                ndm[[x, y]] = self.distance_matrix[[*ix, *iy]];
            }
        }

        return TreeDistanceMatrix::new_from_matrix_and_leaves(ndm, leaves.to_vec());
    }

    fn filter_uniprot_ids(identifiers : Vec<String>) -> Vec<String> {
        // Official uniprot identifier regex
        let re = Regex::new(r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}").unwrap();
        identifiers.into_iter().filter(|id| re.is_match(id)).collect()
    }

    fn uniprot_idx(identifiers : Vec<String>) -> Vec<usize> {
        // Official uniprot identifier regex
        let re = Regex::new(r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}").unwrap();
        identifiers.iter().enumerate().filter(
            |id| re.is_match(id.1)
        ).map(| id | id.0).collect()
    }

    /// Returns the uniprot identifiers in this distance matrix.
    pub fn get_uniprot_ids(&self) -> Vec<String> {
        TreeDistanceMatrix::filter_uniprot_ids(self.leaf_map.keys().map(|k| k.to_string()).collect())
    }

    /// Returns the indices of the uniprot identifiers in this distance matrix.
    pub fn get_uniprot_idx(&self) -> Vec<usize> {
        TreeDistanceMatrix::uniprot_idx(self.leaf_map.keys().map(|k| k.to_string()).collect())
    }

    /// Creates a new distance matrix instance based on parts from the Tree module.
    pub fn new(
        leaf_distance_matrix: Array2<f32>,
        identity_matrix: Array2<usize>,
        leaf_map: HashMap<String, usize>,
    ) -> TreeDistanceMatrix {
        let n_leaves = leaf_map.len();

        let mut leaf_map_inv: HashMap<usize, String> = HashMap::new();
        for (key, value) in leaf_map.iter() {
             leaf_map_inv.insert(*value, key.to_string());
        }

        let distance_matrix = TreeDistanceMatrix::generate_full_distance_matrix(
            leaf_distance_matrix,
            identity_matrix,
            n_leaves,
        );

        TreeDistanceMatrix {
            leaf_map,
            leaf_map_inv,
            distance_matrix,
        }
    }

    /// Performs merging of organisms based on uniprot IDs.
    pub fn merge_organisms(&self, tax_id_map : &BTreeMap<[u8; 10], u32>) -> TreeDistanceMatrix {
        /* Means the distances of all the organisms. */
        let uniprot_ids = self.get_uniprot_ids();
        let uniprot_idx = self.get_uniprot_idx();

        let input_matrix = self.distance_matrix.select(Axis(0), uniprot_idx.as_slice());
        let input_matrix = input_matrix.select(Axis(1), uniprot_idx.as_slice());

        let mut tax_ids = batch_id_tax_mapping(uniprot_ids, tax_id_map);
        let all_tax_ids = tax_ids.clone();

        tax_ids.sort_unstable();
        tax_ids.dedup(); // Deduplicate

        // Allocate a new matrix with the new organisms required
        let mut col_mat : Array2<f32> = Array2::zeros((tax_ids.len(), all_tax_ids.len()));
        let mut final_mat : Array2<f32> = Array2::zeros((tax_ids.len(), tax_ids.len()));

        let enum_tax_ids : Vec<(usize, u32)>= tax_ids.into_iter().enumerate().collect();

        // Mean all the rows in the distance matrix
        for (i, tax_id) in &enum_tax_ids {
            let idx = utils::get_indices_vec(&all_tax_ids, *tax_id);
            let row = input_matrix.select(
                Axis(0),
                idx.as_slice()
            ).mean_axis(Axis(0));

            col_mat.slice_mut(s![*i, ..]).assign(&row);
        }
        
        // Mean all the columns in the distance matrix.
        for (i, tax_id) in &enum_tax_ids {
            let idx = utils::get_indices_vec(&all_tax_ids, *tax_id);
            let row = col_mat.select(Axis(1), idx.as_slice()).mean_axis(Axis(1));
            final_mat.slice_mut(s![*i, ..]).assign(&row);
        }

        return TreeDistanceMatrix::new_from_matrix_and_leaves(
            final_mat,
            enum_tax_ids.iter().map(| i | i.1.to_string()).collect()
        );

    }

    /// Creates a new TreeDistanceMatrix instance based only on a distance matrix and a list of leaves.
    pub fn new_from_matrix_and_leaves(
        distance_matrix: Array2<f32>,
        leaves: Vec<String>,
    ) -> TreeDistanceMatrix {
        let mut leaf_map: HashMap<String, usize> = HashMap::new();
        let mut leaf_map_inv: HashMap<usize, String> = HashMap::new();

        for (i, leaf) in leaves.iter().enumerate() {
            leaf_map.insert(leaf.clone(), i);
            leaf_map_inv.insert(i, leaf.clone());
        }

        return TreeDistanceMatrix {
            leaf_map,
            leaf_map_inv,
            distance_matrix,
        };
    }
}
