use ndarray::prelude::*;
use std::collections::HashMap;
use std::f32;
use tree_distance_matrix::*;

#[derive(Debug)]
pub struct Tree {
    pub leaves: Vec<String>,
    pub branches: Vec<Tree>,
    pub leaf_distances: Vec<f32>,
    pub branch_distances: Vec<f32>,
    levels_from_root: usize,
    branch_number: usize,
}

// Implement partial equality. Trees can have the same structure but be defined differently.
// Thus this equality is intended to be a close approximation. If this returns true, the trees are
// not necessarily exact equals, but they do contain all the same leaves and all siblings are identical.
// http://evolution.genetics.washington.edu/phylip/newicktree.html
impl PartialEq for Tree {
    fn eq(&self, other: &Tree) -> bool {
        let mut leaves_self: Vec<Vec<String>> = Vec::new();
        let mut leaves_other: Vec<Vec<String>> = Vec::new();

        for branch in self.traverse_children() {
            let mut l = branch.0.leaves.clone();
            l.sort();
            leaves_self.push(l);
        }
        for branch in other.traverse_children() {
            let mut l = branch.0.leaves.clone();
            l.sort();
            leaves_other.push(l);
        }

        leaves_self.sort();
        leaves_other.sort();

        return leaves_self == leaves_other;
    }
}

impl Clone for Tree {
    fn clone(&self) -> Tree {
        Tree {
            leaves: self.leaves.iter().map( | l | l.to_string() ).collect(),
            branches: self.branches.iter().map( | b | b.clone() ).collect(),
            leaf_distances: self.leaf_distances.iter().map( | ld | *ld ).collect(),
            branch_distances: self.branch_distances.iter().map( | bd | *bd ).collect(),
            levels_from_root: self.levels_from_root,
            branch_number: self.branch_number
        }
    }
}

impl Tree {
    #[allow(dead_code)]
    pub fn get_leaves(&self) -> HashMap<&String, &f32> {
        /* Returns a HashMap of leaves on this tree with their respective distances. */
        let mut leaves_map = HashMap::new();
        for l in 0..self.leaves.len() {
            leaves_map.insert(
                self.leaves.get(l).expect("Leaves state is invalid"),
                self.leaf_distances
                    .get(l)
                    .expect("Leaves state is invalid: Leaf distance missing."),
            );
        }
        leaves_map
    }

    #[allow(dead_code)]
    pub fn add_leaf(&mut self, name: String, distance: f32) {
        self.leaves.push(name);
        self.leaf_distances.push(distance);
    }

    #[allow(dead_code)]
    pub fn add_branch(&mut self, mut branch: Tree, distance: f32) {
        branch.levels_from_root = self.levels_from_root + 1;
        branch.branch_number = self.branches.len() + 1;

        self.branches.push(branch);
        self.branch_distances.push(distance);
    }

    #[allow(dead_code)]
    fn build_depth_first_path(&self) -> Vec<(&Tree, &f32)> {
        let mut path = Vec::<(&Tree, &f32)>::new();
        for b in 0..self.branches.len() {
            let branch = self.branches.get(b).expect("Branch state invalid");
            let branch_distance = self.branch_distances.get(b).expect(&format!(
                "Branch state is invalid: Branch distance for branch {} missing",
                b
            ));

            path.push((branch, branch_distance));
            if branch.branches.len() > 0 {
                let mut subpaths: Vec<(&Tree, &f32)> = branch.build_depth_first_path();
                path.append(&mut subpaths);
            }
        }
        path
    }

    #[allow(dead_code)]
    pub fn traverse_children(&self) -> Vec<(&Tree, &f32)> {
        let children = self.build_depth_first_path();
        children
    }

    #[allow(dead_code)]
    pub fn traverse_children_snapshot(&self) -> Vec<(Tree, f32)> {
        let children = self.build_depth_first_path();
        let snapshot_children = children.iter().map(| x | (x.0.clone(), *x.1) ).collect();
        snapshot_children
    }
    #[allow(dead_code)]
    fn parse_tree_from_string(
        tree_string: String,
        depth: usize,
        branch_number: usize,
    ) -> (Tree, usize) {
        /* Recursive method that parses a tree structure from a Newick formatted tree. */
        let mut leaves = Vec::<String>::new();
        let mut branches = Vec::<Tree>::new();
        let mut leaf_distances = Vec::<f32>::new();
        let mut branch_distances = Vec::<f32>::new();

        let mut current_leaf = String::new();
        let mut current_distance = String::new();
        let mut read_mode = "LEAF"; // The read mode can be LEAF, DIST, BDONE

        let mut c = 0;
        let mut chr: char;

        while c < tree_string.len() {
//            print!("{}", c);
//            for _c in c.to_string().chars() {
//                print!("{}", (8u8 as char));
//            }
//            print!("{}", (8u8 as char));

            chr = tree_string.chars().nth(c).unwrap();
            c += 1;
//            if depth < 10 {
//                print!("\r");
//                print!("{} chars left.", tree_string.len());
//            }

            if chr == '(' {
                let remainder = tree_string.chars().skip(c).collect();
                let (branch, new_skip) =
                    Tree::parse_tree_from_string(remainder, depth + 1, branches.len());
                c += new_skip;
                branches.push(branch);
                branch_distances.push(f32::NAN);
                read_mode = "BDON";
                continue;
            } else if chr == ')' {
                if read_mode == "DIST" {
                    if current_leaf.len() > 0 {
                        // We found a leaf with a distance.
                        leaves.push(current_leaf.clone().trim().to_string());
                        leaf_distances.push(current_distance.parse().unwrap());
                    } else {
                        // We found the distance of a branch, replace NAN with true distance
                        let bdsize = branch_distances.len();
                        branch_distances[bdsize - 1] = current_distance.parse::<f32>().unwrap();
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
                        branch_distances[bdsize - 1] = current_distance.parse::<f32>().unwrap();
                    }
                } else if read_mode == "LEAF" {
                    // We found a leaf without a distance.
                    leaves.push(current_leaf.clone());
                } /*else if read_mode == "BDON" {
                      read_mode = "LEAF";
                  }*/

                read_mode = "LEAF";
                current_leaf = String::new();
                current_distance = String::new();
            } else {
                if read_mode == "LEAF" {
                    current_leaf.push(chr);
                } else if read_mode == "DIST" {
                    current_distance.push(chr);
                } else if read_mode == "BDON" {
                    continue; // Ignore edge values.
                }
            }
        }

        let tree = Tree {
            leaves: leaves,
            branches: branches,
            leaf_distances: leaf_distances,
            branch_distances: branch_distances,
            levels_from_root: depth,
            branch_number: branch_number,
        };

        return (tree, c);
    }

    #[allow(dead_code)]
    pub fn parse(tree_string: String) -> Tree {
        return Tree::parse_tree_from_string(tree_string, 0, 0).0;
    }

    #[allow(dead_code)]
    pub fn add_root_levels(&mut self, depth: usize) {
        self.levels_from_root = depth;
        for mut branch in self.branches.iter_mut() {
            branch.add_root_levels(depth + 1);
        }
    }

    #[allow(dead_code)]
    pub fn new(
        leaves: Vec<String>,
        branches: Vec<Tree>,
        leaf_distances: Vec<f32>,
        branch_distances: Vec<f32>,
    ) -> Tree {
        if (leaf_distances.len() > 0) & (leaf_distances.len() != leaves.len()) {
            panic!("Leaf distances must be empty or correspond to the number of leaves.");
        }

        if (branch_distances.len() > 0) & (branch_distances.len() != branches.len()) {
            panic!("Branch distances must be empty or correspond to the number of branches.");
        }

        if leaf_distances.len() == 0 {
            let mut leaf_distances: Vec<f32> = Vec::new();
            for _n in 0..branches.len() {
                leaf_distances.push(f32::NAN);
            }
        }

        if branch_distances.len() == 0 {
            let mut branch_distances: Vec<f32> = Vec::new();
            for _n in 0..branches.len() {
                branch_distances.push(f32::NAN);
            }
        }

        return Tree {
            leaves,
            branches,
            leaf_distances,
            branch_distances,
            levels_from_root: 0,
            branch_number: 0,
        };
    }

    #[allow(dead_code)]
    pub fn to_distance_matrix(&self) -> TreeDistanceMatrix {
        let mut n_leaves: usize = 0;
        let mut max_depth = 0;
        let mut leaf_map: HashMap<String, usize> = HashMap::new();

        for child in self.traverse_children() {
            for leaf in child.0.leaves.iter() {
                leaf_map.insert(leaf.clone(), n_leaves);
                n_leaves += 1;
            }
            if child.0.levels_from_root > max_depth {
                max_depth = child.0.levels_from_root;
            }
        }

        let mut leaf_distance_matrix: Array2<f32> = Array2::zeros((n_leaves, max_depth + 1));
        let mut identity_matrix: Array2<usize> = Array2::zeros((n_leaves, max_depth + 1));

        let mut finished_leaves: Vec<usize> = Vec::new();
        let mut accumulated_distances: Vec<f32> = Vec::new();
        let mut internal_nodes: Vec<usize> = Vec::new();

        let mut previous_level = 0;
        let mut previous_branch_number = 0;

        let mut current_level;
        let mut current_branch_number;

        for (c, child) in self.traverse_children().iter().enumerate() {
            let current_tree = child.0;
            let branch_distance = child.1;

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

            for (i, leaf) in current_tree.leaves.iter().enumerate() {
                let leaf_id = leaf_map[leaf];
                for l in 0..current_level {
                    leaf_distance_matrix[[leaf_id, l]] =
                        accumulated_distances[0..l + 1].iter().sum();
                    identity_matrix[[leaf_id, l]] = internal_nodes[l];
                }
                // Add the final accumulated leaf distance
                leaf_distance_matrix[[leaf_id, current_level]] = leaf_distance_matrix
                    [[leaf_id, current_level - 1]]
                    + current_tree.leaf_distances[i];
                finished_leaves.push(leaf_id);
            }

            previous_level = current_level;
            previous_branch_number = current_branch_number;
        }

        // println!("{:?}", identity_matrix);
        // println!("{:?}", leaf_distance_matrix);

        let distance_matrix =
            TreeDistanceMatrix::new(leaf_distance_matrix, identity_matrix, leaf_map);
        distance_matrix
    }
}
