use std::f64;
use std::collections::HashMap;
use ndarray::prelude::*;
use std::process;

#[derive(Debug)]
pub struct TreeDistanceMatrix {
    leaf_map: HashMap<String, usize>,
    distance_matrix: Array2<f64>
}

impl TreeDistanceMatrix {
    #[allow(dead_code)]
    pub fn get_distance(&self, leaf_one : String, leaf_two : String) -> f64 {
        let l1 : usize = self.leaf_map[&leaf_one];
        let l2 : usize = self.leaf_map[&leaf_two];

        self.distance_matrix[[l1, l2]]
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
        let tdm = TreeDistanceMatrix {
            leaf_map,
            distance_matrix : TreeDistanceMatrix::generate_full_distance_matrix(
                leaf_distance_matrix,
                identity_matrix,
                n_leaves
            )
        };
        tdm
    }
}


#[derive(Debug)]
pub struct Tree {
    pub leaves: Vec<String>,
    pub branches: Vec<Tree>,
    pub leaf_distances: Vec<f64>,
    pub branch_distances: Vec<f64>,
    levels_from_root: usize,
    branch_number: usize

}

impl Tree {
    #[allow(dead_code)]
    pub fn get_leaves(&self) -> HashMap<&String, &f64> {
        let mut leaves_map = HashMap::new();
        for l in 0..self.leaves.len() {
            leaves_map.insert(
                self.leaves.get(l).expect("Leaves state is invalid"),
                self.leaf_distances.get(l).expect(
                    "Leaves state is invalid: Leaf distance missing."
                )
            );
        }
        leaves_map
    }

    #[allow(dead_code)]
    pub fn add_leaf(&mut self, name: String, distance: f64) {
        self.leaves.push(name);
        self.leaf_distances.push(distance);
    }

    #[allow(dead_code)]
    pub fn add_branch(&mut self, mut branch: Tree, distance: f64) {
        branch.levels_from_root = self.levels_from_root + 1;
        branch.branch_number = self.branches.len() + 1;

        self.branches.push(branch);
        self.branch_distances.push(distance);
    }

    #[allow(dead_code)]
    fn build_depth_first_path(&self) -> Vec<Vec<(&Tree, &f64)>> {
        let mut path = Vec::<Vec<(&Tree, &f64)>>::new();
        for b in 0..self.branches.len() {
            let branch = self.branches.get(b).expect("Branch state invalid");
            let branch_distance = self.branch_distances.get(b).expect(
                &format!("Branch state is invalid: Branch distance for branch {} missing", b)
            );
            let mut branch_map= Vec::<(&Tree, &f64)>::new();
            branch_map.push((branch, branch_distance));

            path.push(branch_map);
            if branch.branches.len() > 0 {
                let mut subpaths : Vec<Vec<(&Tree, &f64)>> = branch.build_depth_first_path();
                path.append(&mut subpaths);
            }
        }
        path
    }

    #[allow(dead_code)]
    pub fn traverse_children(&self) -> Vec<Vec<(&Tree, &f64)>> {
        let children = self.build_depth_first_path();
        children
    }

    fn parse_tree_from_string(tree_string : String, depth : usize, branch_number : usize) -> (Tree, usize) {
        /* Recursive method that parses a tree structure from a Newick formatted tree. */
        let mut leaves = Vec::<String>::new();
        let mut branches = Vec::<Tree>::new();
        let mut leaf_distances = Vec::<f64>::new();
        let mut branch_distances = Vec::<f64>::new();

        let mut current_leaf = String::new();
        let mut current_distance = String::new();
        let mut read_mode = "LEAF";

        let mut c = 0;
        let mut chr : char;

        while c < tree_string.len() {
            chr = tree_string.chars().nth(c).unwrap();
            c += 1;

            if chr == '(' {
                let remainder = tree_string.chars().skip(c ).collect();
                let (branch, new_skip) = Tree::parse_tree_from_string(
                    remainder,
                    depth + 1,
                    branches.len()
                );
                c += new_skip;
                branches.push(branch);
                branch_distances.push(f64::NAN);
                read_mode = "LEAF";
                continue;

            } else if chr == ')' {
                if read_mode == "DIST" {
                    if current_leaf.len() > 0 {
                        // We found a leaf with a distance.
                        leaves.push(current_leaf.clone().trim().to_string());
                        leaf_distances.push(current_distance.parse().unwrap());
                    } else {
                        // We found the distance of a branch, eplace NAN with true distance
                        let bdsize = branch_distances.len();
                        branch_distances[bdsize - 1] = current_distance.parse::<f64>().unwrap();
                    }
                } else if read_mode == "LEAF" {
                    // We found a leaf without a distance.
                    leaves.push(current_leaf.clone().trim().to_string());
                }
                break;

            } else if chr == ':' {
                read_mode = "DIST";
            } else if chr == ',' {
                // End of an item in a list
                if read_mode == "DIST" {
                    if current_leaf.len() > 0 {
                        // We found a leaf with a distance.
                        leaves.push(current_leaf.clone());
                        leaf_distances.push(current_distance.parse().unwrap());
                    } else {
                        // We found the distance of a branch, replace NAN with true distance
                        let bdsize = branch_distances.len();
                        branch_distances[bdsize - 1] = current_distance.parse::<f64>().unwrap();
                    }
                } else if read_mode == "LEAF" {
                    // We found a leaf without a distance.
                    leaves.push(current_leaf.clone());
                }

                read_mode = "LEAF";
                current_leaf = String::new();
                current_distance = String::new();

            } else {
                if read_mode == "LEAF" {
                    current_leaf.push(chr);
                } else if read_mode == "DIST" {
                    current_distance.push(chr);
                }
            }
        }

        let tree = Tree {
            leaves: leaves,
            branches: branches,
            leaf_distances: leaf_distances,
            branch_distances: branch_distances,
            levels_from_root: depth,
            branch_number: branch_number
        };

        return (tree, c);
    }

    pub fn parse(tree_string : String) -> Tree {
        return Tree::parse_tree_from_string(tree_string, 0, 0).0;
    }

    pub fn new(leaves : Vec<String>, branches : Vec<Tree>, leaf_distances : Vec<f64>, branch_distances : Vec<f64>) -> Tree {

        if (leaf_distances.len() > 0) & (leaf_distances.len() != leaves.len()) {
            panic!("Leaf distances must be empty or correspond to the number of leaves.");
        }

        if (branch_distances.len() > 0) & (branch_distances.len() != branches.len()) {
            panic!("Branch distances must be empty or correspond to the number of branches.");
        }

        if leaf_distances.len() == 0 {
            let mut leaf_distances : Vec<f64> = Vec::new();
            for _n in 0..branches.len() {
                leaf_distances.push(f64::NAN);
            }
        }

        if branch_distances.len() == 0 {
            let mut branch_distances : Vec<f64> = Vec::new();
            for _n in 0..branches.len() {
                branch_distances.push(f64::NAN);
            }
        }

        return Tree {
            leaves,
            branches,
            leaf_distances,
            branch_distances,
            levels_from_root: 0,
            branch_number: 0
        };
    }

    pub fn to_distance_matrix(&self)-> TreeDistanceMatrix {
        let mut n_leaves : usize = 0;
        let mut max_depth = 0;
        let mut leaf_map : HashMap<String, usize> = HashMap::new();

        for child in self.traverse_children() {
            for leaf in child[0].0.leaves.iter() {
                leaf_map.insert(leaf.clone(), n_leaves);
                n_leaves += 1;
            }
            if child[0].0.levels_from_root > max_depth {
                max_depth = child[0].0.levels_from_root;
            }
        }

        let mut leaf_distance_matrix : Array2<f64> = Array2::zeros((n_leaves, max_depth + 1));
        let mut identity_matrix : Array2<usize> = Array2::zeros((n_leaves, max_depth + 1));

        let mut finished_leaves : Vec<usize> = Vec::new();
        let mut accumulated_distances : Vec<f64> = Vec::new();
        let mut internal_nodes : Vec<usize> = Vec::new();

        let mut previous_level = 0;
        let mut previous_branch_number = 0;

        let mut current_level;
        let mut current_branch_number;

        for (c, child) in self.traverse_children().iter().enumerate() {
            let current_tree = child[0].0;
            let branch_distance = child[0].1;

            current_level = current_tree.levels_from_root;
            current_branch_number = current_tree.branch_number;

            if previous_level < current_level {
                for _l in current_level..previous_level {
                    accumulated_distances.pop();
                    internal_nodes.pop();
                }
            } else if previous_branch_number < current_branch_number {
                accumulated_distances.pop();
                internal_nodes.pop();
            }

            if (*branch_distance).is_nan() {
                accumulated_distances.push(0.0);
            } else {
                accumulated_distances.push(*branch_distance);
            }

            internal_nodes.push(c);

            for (i, leaf)in current_tree.leaves.iter().enumerate() {
                let leaf_id = leaf_map[leaf];
                for l in 0..current_level {
                    leaf_distance_matrix[[leaf_id, l]] = accumulated_distances[0..l + 1].iter().sum();
                    identity_matrix[[leaf_id, l]] = internal_nodes[l];
                }
                // Add the final accumulated leaf distance
                leaf_distance_matrix[[leaf_id, current_level]] =
                    leaf_distance_matrix[[leaf_id, current_level - 1]] + current_tree.leaf_distances[i];
                finished_leaves.push(leaf_id);
            }

            previous_level = current_level;
            previous_branch_number = current_branch_number;
        }

        println!("{:?}", identity_matrix);
        println!("{:?}", leaf_distance_matrix);

        let distance_matrix = TreeDistanceMatrix::new(
            leaf_distance_matrix,
            identity_matrix,
            leaf_map
        );
        distance_matrix
    }
}


