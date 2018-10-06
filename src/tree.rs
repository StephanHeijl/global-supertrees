use std::f64;
use std::collections::HashMap;
use ndarray::prelude::*;


#[derive(Debug)]
pub struct Tree {
    pub leaves: Vec<String>,
    pub branches: Vec<Tree>,
    pub leaf_distances: Vec<f64>,
    pub branch_distances: Vec<f64>
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
    pub fn add_branch(&mut self, branch: Tree, distance: f64) {
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

    fn parse_tree_from_string(tree_string : String) -> (Tree, usize) {
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
                let (branch, new_skip) = Tree::parse_tree_from_string(remainder);
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
            branch_distances: branch_distances
        };

        return (tree, c);
    }

    pub fn parse(tree_string : String) -> Tree {
        return Tree::parse_tree_from_string(tree_string).0;
    }

    pub fn to_distance_matrix(&self)-> Array2<f64> {
        let mut n_leaves = 0;
        let mut n_non_leaf_nodes = 0;
        for child in self.traverse_children() {
            //println!("{:?}", child);
            n_non_leaf_nodes += 1;
            n_leaves += child[0].0.leaves.len();
        }

        let distance_matrix = Array::zeros((n_leaves, n_leaves).f());


        distance_matrix
    }
}


