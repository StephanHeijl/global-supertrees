use std::f64;
use std::collections::HashMap;
use ndarray::prelude::*;
use tree::*;

#[derive(Debug)]
pub struct TreeDistanceMatrix {
    pub leaf_map: HashMap<String, usize>,
    pub leaf_map_inv: HashMap<usize, String>,
    pub distance_matrix: Array2<f64>
}

impl TreeDistanceMatrix {
    pub fn calculate_q_matrix(dm : &Array2<f64>) -> Array2<f64> {
        let mut q_matrix : Array2<f64 >= Array2::zeros(dm.raw_dim());
        let n = dm.shape()[0];  // Number of taxa

        for i in 0..n {
            for j in 0..n {
                q_matrix[[i, j]] =  ((n as f64 - 2.0) * dm[[i, j]]) -
                                     dm.row(i).scalar_sum() -
                                     dm.column(j).scalar_sum();
            }
        }
        q_matrix
    }

    pub fn calculate_pair_distance_matrix(i : usize, j : usize, dm : &Array2<f64>) -> Array2<f64> {
        // Create new distance matrix with the new pair as a single node.
        let mut new_distances : Array2<f64 >= Array2::zeros(
            (dm.shape()[0] - 1, dm.shape()[1] - 1)
        );
        // New indexes
        let mut x_n = 0;
        let mut y_n;
        let final_idx = new_distances.shape()[0] - 1;

        for x in 0..dm.shape()[0] {
            y_n = 0;
            for y in 0..dm.shape()[1] {
                //println!("{},{}, {:?}, {:?}, {}, {}, {}, {}", x, y, new_distances.shape(), dm.shape(), x_n, y_n, i, j);
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

    pub fn calculate_new_leaf_distances(f : usize, g : usize, dm : &Array2<f64>) -> (f64, f64) {
        /* Finds the leaf distances for the pair inside the Tree node. */
        // f and g are the paired taxa and u is the new created node.
        let n = dm.shape()[0] as f64;

        let d_fu = 0.5 * dm[[f, g]] + (1.0 / (2.0 * (n - 2.0))) * (dm.row(f).scalar_sum() - dm.row(g).scalar_sum());
        let d_gu = dm[[f, g]] - d_fu;

        return (d_fu, d_gu);
    }

    pub fn find_min_matrix_ne(matrix : Array2<f64>) -> [usize; 2] {
        /* Returns the lowest index in the matrix where i != j */
        let mut lowest : f64 = f64::MAX;
        let mut idx = [0, 0];
        let mut v : f64;
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

        let mut trees : Vec<Tree> = Vec::new();
        let mut iteration = 0;

        while distance_matrix.shape()[0] > 2 {
            let q_matrix = TreeDistanceMatrix::calculate_q_matrix(
                &distance_matrix
            );
            let lowest_pair = TreeDistanceMatrix::find_min_matrix_ne(q_matrix);

            let pair_distance_matrix = TreeDistanceMatrix::calculate_pair_distance_matrix(
                lowest_pair[0], lowest_pair[1], &distance_matrix
            );
            let pair_leaf_distances = TreeDistanceMatrix::calculate_new_leaf_distances(
                lowest_pair[0], lowest_pair[1], &distance_matrix
            );

            let tree_idx : Vec<usize> = ((pair_distance_matrix.shape()[0] - iteration)..(pair_distance_matrix.shape()[0] + 1)).collect();
            println!("{:?}", tree_idx);

            if tree_idx.contains(&lowest_pair[0]) {

            }

            let mut leaves = Vec::new();
            leaves.push(self.leaf_map_inv[&lowest_pair[0]].clone());
            leaves.push(self.leaf_map_inv[&lowest_pair[1]].clone());

            let inner_tree = Tree::new(
                leaves,
                vec!(),
                vec!(pair_leaf_distances.0, pair_leaf_distances.1),
                vec!()
            );
            distance_matrix = pair_distance_matrix;



            trees.push(inner_tree);

            let new_pair_dest =

            println!("{:?} -> {:?} ({})", distance_matrix.shape(), lowest_pair, iteration);
            iteration += 1;
        }

        return Tree::new(vec!(), vec!(), vec!(), vec!());
    }


    #[allow(dead_code)]
    pub fn get_distance(&self, leaf_one : String, leaf_two : String) -> f64 {
        let l1 : usize = self.leaf_map[&leaf_one];
        let l2 : usize = self.leaf_map[&leaf_two];

        self.distance_matrix[[l1, l2]]
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


    fn generate_partial_distance_matrix(leaf_distance_matrix : &Array2<f64>,
                                        identity_matrix: &Array2<usize>,
                                        start: usize,
                                        stop: usize) -> Array2<f64> {

        let n : usize;
        if stop >= leaf_distance_matrix.shape()[0] {
            n = leaf_distance_matrix.shape()[0] - 1 - start;
        } else {
            n = stop - start;
        }

        //println!("{}->{}  {:?}", start, stop, leaf_distance_matrix.shape());
        let mut distance_matrix : Array2<f64 >= Array2::zeros((n, n));

        for x in 0..n {
            for y in 0..n {
                if x == y {
                    distance_matrix[[x, y]] = 0.0;
                    continue;
                }
                let x_full = x + start;
                let y_full = y + start;

                let (xi, yi) = TreeDistanceMatrix::find_final_parent(x_full, y_full, identity_matrix);
                let fca = TreeDistanceMatrix::find_first_common_ancestor(x_full, y_full, identity_matrix);

                // Find the total distance from each leaf to the root.
                let leaf_root_distance : f64 = leaf_distance_matrix[[x_full, xi]] + leaf_distance_matrix[[y_full, yi]];

                // Find the distance between the root and the first common ancestor.
                let root_fca_distance: f64  = leaf_distance_matrix[[x_full, fca]];
                let distance = f64::abs(leaf_root_distance - (root_fca_distance * 2.0));
                distance_matrix[[x, y]] = distance;
            }
        }

        distance_matrix
    }


    fn generate_full_distance_matrix(leaf_distance_matrix : Array2<f64>,
                                     identity_matrix: Array2<usize>,
                                     n_leaves : usize) -> Array2<f64> {
        /* Generates a full distance matrix */

        //let mut distance_matrices : Vec<Array2<f64>> = Vec::new();
        let max_size = 32.0;
        let parts = ((n_leaves as f64) / max_size).ceil() as usize;

        let starts : Vec<usize> = (0..parts).collect();
        println!("{:?}", starts);

        let distance_matrices : Vec<Array2<f64>> = starts.iter().map(
            |&s| TreeDistanceMatrix::generate_partial_distance_matrix(
                &leaf_distance_matrix,
                &identity_matrix,
                (s * max_size as usize),
                ((s+1) * max_size as usize)
            )
        ).collect();

        println!("{:?}", distance_matrices);

        let mut distance_matrix : Array2<f64 >= Array2::zeros((2, 2));
        distance_matrix
    }

    pub fn new(leaf_distance_matrix : Array2<f64>,
               identity_matrix: Array2<usize>,
               leaf_map: HashMap<String, usize>) -> TreeDistanceMatrix {
        let n_leaves = leaf_map.len();

        let mut leaf_map_inv : HashMap<usize, String> = HashMap::new();
        for (key, value) in leaf_map.iter() {
            leaf_map_inv.insert(*value, key.to_string());
        }

        let tdm = TreeDistanceMatrix {
            leaf_map,
            leaf_map_inv,
            distance_matrix : TreeDistanceMatrix::generate_full_distance_matrix(
                leaf_distance_matrix,
                identity_matrix,
                n_leaves
            )
        };
        tdm
    }
}