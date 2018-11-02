use std::f64;
use std::collections::HashMap;
use ndarray::prelude::*;

#[derive(Debug)]
pub struct TreeDistanceMatrix {
    pub leaf_map: HashMap<String, usize>,
    pub leaf_map_inv: HashMap<usize, String>,
    pub distance_matrix: Array2<f64>
}

impl TreeDistanceMatrix {
    #[allow(dead_code)]
    pub fn get_distance(&self, leaf_one : String, leaf_two : String) -> f64 {
        let l1 : usize = self.leaf_map[&leaf_one];
        let l2 : usize = self.leaf_map[&leaf_two];

        self.distance_matrix[[l1, l2]]
    }

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

    pub fn shape(&self) -> &[usize]{
        return self.distance_matrix.shape();
    }

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

    fn generate_full_distance_matrix(leaf_distance_matrix : Array2<f64>,
                                     identity_matrix: Array2<usize>,
                                     n_leaves : usize) -> Array2<f64> {
        /* Generates a full distance matrix */

        let mut distance_matrix : Array2<f64 >= Array2::zeros((n_leaves, n_leaves));

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
                let leaf_root_distance : f64 = leaf_distance_matrix[[x, xi]] + leaf_distance_matrix[[y, yi]];

                // Find the distance between the root and the first common ancestor.
                let root_fca_distance: f64  = leaf_distance_matrix[[x, fca]];
                let distance = f64::abs(leaf_root_distance - (root_fca_distance * 2.0));
                distance_matrix[[x, y]] = distance;
            }
        }

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