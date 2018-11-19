use std::f32;
use std::collections::HashMap;
use ndarray::prelude::*;
use tree::*;
use ndarray::stack;
use std::cmp::min;
use rayon::prelude::*;
use utils;


#[derive(Debug)]
pub struct TreeDistanceMatrix {
    pub leaf_map: HashMap<String, usize>,
    pub leaf_map_inv: HashMap<usize, String>,
    pub distance_matrix: Array2<f32>
}

impl TreeDistanceMatrix {
    pub fn calculate_q_matrix(dm : &Array2<f32>) -> Array2<f32> {
        let mut q_matrix : Array2<f32 >= Array2::zeros(dm.raw_dim());
        let n = dm.shape()[0];  // Number of taxa

        for i in 0..n {
            for j in 0..n {
                q_matrix[[i, j]] =  ((n as f32 - 2.0) * dm[[i, j]]) -
                    dm.row(i).scalar_sum() -
                    dm.column(j).scalar_sum();
            }
        }
        q_matrix
    }

    pub fn calculate_pair_distance_matrix(i : usize, j : usize, dm : &Array2<f32>) -> Array2<f32> {
        // Create new distance matrix with the new pair as a single node.
        let mut new_distances : Array2<f32 >= Array2::zeros(
            (dm.shape()[0] - 1, dm.shape()[1] - 1)
        );
        // New indexes
        let mut x_n = 0;
        let mut y_n;
        let final_idx = new_distances.shape()[0] - 1;

        for x in 0..dm.shape()[0] {
            y_n = 0;
            for y in 0..dm.shape()[1] {
                if (x == i) | (x == j) {
                    new_distances[[final_idx, y_n]] = 0.5 * (
                        dm[[i, y]] + dm[[j, y]] - dm[[i, j]]
                    );
                    if (y != i) & (y != j) {
                        y_n += 1;
                    }
                } else if (y == i) | (y == j) {
                    new_distances[[x_n, final_idx]] = 0.5 * (
                        dm[[i, x]] + dm[[j, x]] - dm[[i, j]]
                    );
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

    pub fn calculate_new_leaf_distances(f : usize, g : usize, dm : &Array2<f32>) -> (f32, f32) {
        /* Finds the leaf distances for the pair inside the Tree node. */

        // f and g are the paired taxa and u is the new created node.
        let n = dm.shape()[0] as f32;

        let d_fu = 0.5 * dm[[f, g]] + (1.0 / (2.0 * (n - 2.0))) * (dm.row(f).scalar_sum() - dm.row(g).scalar_sum());
        if d_fu.is_nan() {
            return (0.0, 0.0);
        }
        let d_gu = dm[[f, g]] - d_fu;

        return (d_fu, d_gu);
    }

    pub fn find_min_matrix_ne(matrix : Array2<f32>) -> [usize; 2] {
        /* Returns the lowest index in the matrix where i != j */
        let mut lowest : f32 = f32::MAX;
        let mut idx = [0, 0];
        let mut v : f32;
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

    pub fn neighbour_joining(&self) -> Tree {
        let mut distance_matrix = self.distance_matrix.clone();
        let mut leaf_map_inv_vec : Vec<(usize, String)> = self.leaf_map_inv.iter().map(
            | val | (*val.0, val.1.clone())
        ).collect();

        leaf_map_inv_vec.sort_by(|a, b| a.0.cmp(&b.0));

        let mut trees : HashMap<usize, Tree> = HashMap::new();


        while distance_matrix.shape()[0] > 1 {
            // Calculate the Q matrix
            let q_matrix = TreeDistanceMatrix::calculate_q_matrix(
                &distance_matrix
            );
            // Find the pair of leaves that have the lowest value in the q_matrix
            let lowest_pair = TreeDistanceMatrix::find_min_matrix_ne(q_matrix);

            // Calculate the distance matrix from each of the taxa to the new node
            // In this pair_distance_matrix, the final row/col is now designated the previously
            // selected pair.
            let pair_distance_matrix = TreeDistanceMatrix::calculate_pair_distance_matrix(
                lowest_pair[0], lowest_pair[1], &distance_matrix
            );

            // Calculate the distances between the new node and the two merged taxa
            let pair_leaf_distances = TreeDistanceMatrix::calculate_new_leaf_distances(
                lowest_pair[0], lowest_pair[1], &distance_matrix
            );

            let mut leaves : Vec<String> = Vec::new();
            let mut branches: Vec<Tree> = Vec::new();
            let mut leaf_distances: Vec<f32> = Vec::new();
            let mut branch_distances: Vec<f32> = Vec::new();

            //println!("{:?}", trees);
            // Check if the first node is already a tree
            let mut rem_leaves : Vec<usize> = Vec::new();
            if trees.contains_key(&lowest_pair[0]) {
                branches.push(trees.remove(&lowest_pair[0]).unwrap());
                branch_distances.push(pair_leaf_distances.0);
            } else {
                rem_leaves.push(lowest_pair[0]);
                leaves.push(leaf_map_inv_vec.get(lowest_pair[0]).unwrap().1.clone());
                leaf_distances.push(pair_leaf_distances.0);
            }

            // Check if the second node is already a tree
            if trees.contains_key(&lowest_pair[1]) {
                branches.push(trees.remove(&lowest_pair[1]).unwrap());
                branch_distances.push(pair_leaf_distances.1);
            } else {
                rem_leaves.push(lowest_pair[1]);
                leaves.push(leaf_map_inv_vec.get(lowest_pair[1]).unwrap().1.clone());
                leaf_distances.push(pair_leaf_distances.1);
            }

            // Remove leaves from the leaf map if they have been selected
            let mut leaf_map_inv_vec_new = Vec::new();
            for (i, el) in leaf_map_inv_vec.clone().iter().enumerate() {
                if !rem_leaves.contains(&i) {
                    leaf_map_inv_vec_new.push(el.clone());
                }
            }
            leaf_map_inv_vec = leaf_map_inv_vec_new.clone();

            // Shift all the trees to conform to the new shifted index, the distance matrix shrunk
            // and the new node is the new final tree.
            let mut trees_new : HashMap<usize, Tree> = HashMap::new();
            // Get the tree keys without borrowing them.
            let tree_keys : Vec<usize> = trees.keys().map(| x | *x).collect();
            for idx in tree_keys {
                trees_new.insert(idx - 2, trees.remove(&idx).unwrap());
            }
            trees = trees_new;

            // Create a new Tree based on the selected pair of nodes
            let inner_tree = Tree::new(
                leaves,
                branches,
                leaf_distances,
                branch_distances
            );
            distance_matrix = pair_distance_matrix;

            trees.insert(distance_matrix.shape()[0] - 1, inner_tree);
        }
        assert_eq!(trees.len(), 1);
        let tree_keys : Vec<usize> = trees.keys().map(| x | *x).collect();
        let mut result_tree = trees.remove(&tree_keys[0]).unwrap();
        result_tree.add_root_levels(0);
        return result_tree;

    }


    #[allow(dead_code)]
    pub fn get_distance(&self, leaf_one : String, leaf_two : String) -> Option<f32> {
        return self.get_distance_ref(&leaf_one, &leaf_two);
    }

    #[allow(dead_code)]
    pub fn get_distance_ref(&self, leaf_one : &String, leaf_two : &String) -> Option<f32>  {
        let l1 : Option<&usize> = self.leaf_map.get(leaf_one);
        let l2 : Option<&usize> = self.leaf_map.get(leaf_two);
        let x : usize;
        let y : usize;
        match l1 {
            Some(i) => {x = *i},
            None => { return None; }
        }
        match l2 {
            Some(i) => {y = *i},
            None => { return None; }
        }
        return Some(self.distance_matrix[[x, y]]);
    }

    #[allow(dead_code)]
    pub fn to_csv(&self) -> String {
        let mut csv = String::from("");
        let delim = ',';

        //println!("{} - {:?}", self.leaf_map_inv.len(), self.distance_matrix.shape());

        let mut max_key: usize = 0;
        for key in self.leaf_map_inv.keys() {
            if key > &max_key {
                max_key = *key;
            }
        }

        let mut leaves : Vec<String> = Vec::new();

        for l in 0..max_key + 1 {
            match self.leaf_map_inv.get(&l) {
                Some(name) => {
                    csv.push_str(name);
                    leaves.push(name.to_string());
                    if l < max_key {
                        csv.push(delim);
                    }
                },
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
                if v < leaves.len() {
                    csv.push(delim);
                }
                v += 1;
            }
            csv.push('\n');
        }
        return csv;
    }

    #[allow(dead_code)]
    pub fn shape(&self) -> &[usize]{
        return self.distance_matrix.shape();
    }

    #[allow(dead_code)]
    pub fn find_first_zero(identity_row : ArrayView1<usize>) -> usize {
        /* Returns the index of the first zero in a 1D ArrayView from position 1.
         Returns 0 if no zero has been found. */
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


    fn find_final_parent(l1 : usize, l2 : usize, identity_matrix : &Array2<usize>) -> (usize, usize) {
        /* Returns the indices of the last zeros for two rows in the identity matrix based on the 2 leaf identifiers (l1, l2). */
        let p1 = TreeDistanceMatrix::find_first_zero(identity_matrix.slice(s![l1, ..]));
        let p2 = TreeDistanceMatrix::find_first_zero(identity_matrix.slice(s![l2, ..]));

        (p1, p2)
    }

    pub fn find_first_common_ancestor(l1 : usize, l2 : usize, identity_matrix : &Array2<usize>) -> usize {
        /* Find the first common ancestor of two leaves.  If they do not share any explicit ancestors,
         the first common ancestor becomes the root of the tree (0). If the leaves are on the same (sub)tree,
         their shared subtree is returned. The ancestor trees are returned as a number representing the subtree id. */

        let id_row_1 = identity_matrix.slice(s![l1, ..]);
        let id_row_2 = identity_matrix.slice(s![l2, ..]);

        if id_row_1 == id_row_2 {
            return TreeDistanceMatrix::find_final_parent(l1, l2, identity_matrix).0 - 1
        }

        for i in 1..id_row_1.len() {
            if id_row_1[i] != id_row_2[i] {
                return i - 1;
            }
        }
        return 0;
    }

    #[cfg(not(feature = "singlecore"))]
    fn generate_partial_distance_matrix(leaf_distance_matrix : &Array2<f32>,
                                        identity_matrix: &Array2<usize>,
                                        iteration : usize,
                                        size : usize ) -> Array2<f32> {



        //println!("{}->{}  {:?}", start, stop, leaf_distance_matrix.shape());
        let mut distance_matrix : Array2<f32> = Array2::zeros((size, size));
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
                let leaf_root_distance : f32 = leaf_distance_matrix[[x, xi]] + leaf_distance_matrix[[y, yi]];

                // Find the distance between the root and the first common ancestor.
                let root_fca_distance: f32  = leaf_distance_matrix[[x, fca]];
                let distance = f32::abs(leaf_root_distance - (root_fca_distance * 2.0));

                // Add the distance to the distance matrix
                distance_matrix[[xd, yd]] = distance;
            }
        }

        distance_matrix
    }

    #[cfg(not(feature = "singlecore"))]
    fn generate_full_distance_matrix(leaf_distance_matrix : Array2<f32>,
                                     identity_matrix: Array2<usize>,
                                     n_leaves : usize) -> Array2<f32> {
        /* Generates a full distance matrix */

        //let mut distance_matrices : Vec<Array2<f32>> = Vec::new();
        let max_size = 2048.0;
        let parts = ((n_leaves as f32) / max_size).ceil() as usize;

        let starts : Vec<usize> = (0..parts.pow(2)).collect();


        let distance_matrices : Vec<Array2<f32>> = starts.par_iter().map(
            |&s| TreeDistanceMatrix::generate_partial_distance_matrix(
                &leaf_distance_matrix,
                &identity_matrix,
                s,
                max_size as usize
            )
        ).collect();

        // Compile partial distance matrices
        let mut rows : Vec<Array2<f32>> = Vec::new();
        for x in 0..parts {
            let start_x = x * parts;
            let stop_x = x * parts + parts;
            let row : Vec<ndarray::ArrayBase<ndarray::ViewRepr<&f32>, ndarray::Dim<[usize; 2]>>> = distance_matrices[start_x..stop_x].iter().map(
                |dm| dm.view()
            ).collect();
            rows.push(stack(Axis(1), &row).unwrap());
        }

        let row_views : Vec<ndarray::ArrayBase<ndarray::ViewRepr<&f32>, ndarray::Dim<[usize; 2]>>> = rows.iter().map(
            |row| row.view()
        ).collect();

        let distance_matrix = stack(Axis(0), &row_views).unwrap();
        distance_matrix.slice(s![..n_leaves, ..n_leaves]).to_owned()
    }

    #[cfg(feature = "singlecore")]
    fn generate_full_distance_matrix(leaf_distance_matrix : Array2<f32>,
                                     identity_matrix: Array2<usize>,
                                     n_leaves : usize) -> Array2<f32> {
        /* Generates a full distance matrix */

        let mut distance_matrix : Array2<f32 >= Array2::zeros((n_leaves, n_leaves));

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
                let leaf_root_distance : f32 = leaf_distance_matrix[[x, xi]] + leaf_distance_matrix[[y, yi]];

                // Find the distance between the root and the first common ancestor.
                let root_fca_distance: f32  = leaf_distance_matrix[[x, fca]];
                let distance = f32::abs(leaf_root_distance - (root_fca_distance * 2.0));
                distance_matrix[[x, y]] = distance;
            }
        }

        distance_matrix
    }

    pub fn get_partial_distance_matrix(&self, leaves : &Vec<String>) -> TreeDistanceMatrix {
        let mut leaf_map : HashMap<String, usize> = HashMap::new();
        let mut leaf_map_inv : HashMap<usize, String> = HashMap::new();
        let self_leaves : Vec<String> = utils::keys_vec(&self.leaf_map);

        let mut l = 0;
        let mut dm_idx : Vec<usize> = Vec::new();
        for leaf in leaves.iter() {
            if self_leaves.contains(&leaf) {
                leaf_map.insert(leaf.clone(), l);
                leaf_map_inv.insert(l, leaf.clone());
                dm_idx.push(self.leaf_map[leaf]);
                l += 1;
            }
        }

        let mut ndm : Array2<f32> = Array2::zeros((l, l));

        for (x, ix) in dm_idx.iter().enumerate() {
            for (y, iy) in dm_idx.iter().enumerate() {
                ndm[[x, y]] = self.distance_matrix[[*ix, *iy]];
            }
        }

        return TreeDistanceMatrix::new_from_matrix_and_leaves(ndm, leaves.to_vec());

    }


    pub fn new(leaf_distance_matrix : Array2<f32>,
               identity_matrix: Array2<usize>,
               leaf_map: HashMap<String, usize>) -> TreeDistanceMatrix {
        let n_leaves = leaf_map.len();

        let mut leaf_map_inv : HashMap<usize, String> = HashMap::new();
        for (key, value) in leaf_map.iter() {
            leaf_map_inv.insert(*value, key.to_string());
        }

        let distance_matrix = TreeDistanceMatrix::generate_full_distance_matrix(
            leaf_distance_matrix,
            identity_matrix,
            n_leaves
        );

        TreeDistanceMatrix {
            leaf_map,
            leaf_map_inv,
            distance_matrix
        }
    }

    pub fn new_from_matrix_and_leaves(distance_matrix : Array2<f32>, leaves : Vec<String>) -> TreeDistanceMatrix {

        let mut leaf_map : HashMap<String, usize> = HashMap::new();
        let mut leaf_map_inv : HashMap<usize, String> = HashMap::new();

        for (i, leaf) in leaves.iter().enumerate() {
            leaf_map.insert(leaf.clone(), i);
            leaf_map_inv.insert(i, leaf.clone());
        }

        return TreeDistanceMatrix {
            leaf_map,
            leaf_map_inv,
            distance_matrix
        };
    }
}