use ndarray::prelude::*;
use std::collections::HashMap;
use std::f32;
use regex::Regex;
use tree_distance_matrix::*;

use petgraph::prelude::NodeIndex;
use petgraph::visit::Dfs;
use petgraph::graph::{node_index, EdgeReference};
use petgraph::{Graph, Incoming};
use petgraph::visit::EdgeRef;
use petgraph::dot::{Dot, Config};

/// A tree structure, based on the pet-graph Graph structure. Using this ensures that there are no
/// loops in the graph and that a tree-like structure is guaranteed.
#[derive(Debug)]
pub struct Tree {
    pub graph : Graph<String, f32>,
}

/// A single level on a tree. Contains leaves with distances and a notion of distance from the root
/// of the tree.
#[derive(Debug)]
pub struct Level {
    pub leaves : Vec<String>,
    pub leaf_nodes : Vec<NodeIndex<u32>>,
    pub level_distance : f32,
    pub leaf_distances : Vec<f32>,
    pub levels_from_root: usize,
}

impl Level {
    pub fn to_newick(&self) -> String {
        let mut newick = String::new();
        for (l, leaf) in self.leaves.iter().enumerate() {
            let dist = self.leaf_distances[l];
            if dist.is_nan() {
                newick.push_str(leaf);
            } else {
                newick.push_str(&format!("{}:{}", leaf, dist));
            }
            if l + 1 < self.leaves.len() {
                newick.push(',');
            }

        }
        return newick;
    }
}

impl PartialEq for Tree {
    fn eq(&self, other: &Tree) -> bool {
        let mut leaves_self: Vec<Vec<String>> = Vec::new();
        let mut leaves_other: Vec<Vec<String>> = Vec::new();

        for level in self.traverse_children() {
            let mut l = level.leaves.clone();
            l.sort();
            leaves_self.push(l);
        }
        for level in other.traverse_children() {
            let mut l = level.leaves.clone();
            l.sort();
            leaves_other.push(l);
        }

        leaves_self.sort();
        leaves_other.sort();

        return leaves_self == leaves_other;
    }
}

/// Functions that relate to pure trees.
impl Tree {
    /// Returns every leaf in the tree on every level.
    pub fn get_leaves(&self) -> Vec<String> {
        let mut leaves : Vec<String> = Vec::new();
        for id in self.graph.node_indices() {
            if let Some(thing) = self.graph.node_weight(id) {
                if !thing.starts_with(">>") {
                    leaves.push(thing.to_string());
                }
            }
        }
        leaves
    }

    fn filter_uniprot_ids(identifiers : Vec<String>) -> Vec<String> {
        // Official uniprot identifier regex
        let re = Regex::new(r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}").unwrap();
        identifiers.into_iter().filter(|id| re.is_match(id)).collect()
    }

    /// Returns the leaves in the tree that conform to the Uniprot identifier standard.
    pub fn get_uniprot_ids(&self) -> Vec<String> {
        Tree::filter_uniprot_ids(self.get_leaves())
    }

    /// Returns the edge that connects this node to its parent.
    pub fn get_parent_edge(&self, node : NodeIndex<u32>) -> Option<EdgeReference<f32, u32>> {
        return self.graph.edges_directed(node, Incoming).nth(0);
    }

    /// Returns the node above this node in the tree structure. (A node that is closer to the root.)
    pub fn get_parent_node(&self, node : NodeIndex<u32>) -> NodeIndex<u32> {
        let parent_edge = self.get_parent_edge(node).expect("Node does not have a parent.");
        parent_edge.source()
    }

    /// Returns a vector of Levels in a depth-first fashion.
    pub fn traverse_children(&self) -> Vec<Level> {

        let mut children : Vec<Level> = Vec::new();
        let mut cl_names : Vec<String> = Vec::new();
        let mut cl_nodes : Vec<NodeIndex<u32>> = Vec::new();
        let mut cl_dist : Vec<f32> = Vec::new();
        let mut previous_branch_distance = f32::NAN;

        let mut dfs = Dfs::new(&self.graph,  node_index(0));

        while let Some(node) = dfs.next(&self.graph) {
            let node_name = self.graph.node_weight(node).unwrap();
            let mut node_distance = f32::NAN;

            match self.get_parent_edge(node) {
                Some(edge) => { node_distance = *edge.weight(); },
                None => { }
            }

            if node_name == "<root>" {
                continue;
            }

            if node_name.starts_with(">>") {
                //println!("Branch node ({}) distance: {}", node_name, node_distance);
                if cl_nodes.len() == 0 {
                    continue;
                }
                let levels_from_root = self.get_levels_from_root(*cl_nodes.get(0).expect(
                    "Found no nodes on current level."
                ));
                let level = Level {
                    leaves : cl_names.clone(),
                    leaf_nodes : cl_nodes.clone(),
                    leaf_distances: cl_dist.clone(),
                    level_distance: previous_branch_distance,
                    levels_from_root: levels_from_root,
                };
                children.push(level);

                cl_names.clear();
                cl_nodes.clear();
                cl_dist.clear();

                // Record the true branch distance at the time the branch was encountered.
                previous_branch_distance = node_distance;
            } else {
                cl_names.push(node_name.to_string());
                cl_nodes.push(node);
                cl_dist.push(node_distance);
            }

            // Add the final level if there are no more nodes in the stack.
            if dfs.stack.len() == 0 {
                let levels_from_root = self.get_levels_from_root(*cl_nodes.get(0).expect(
                    "Found no nodes on current level."
                ));

                let level = Level {
                    leaves: cl_names.clone(),
                    leaf_nodes: cl_nodes.clone(),
                    leaf_distances: cl_dist.clone(),
                    level_distance: previous_branch_distance,
                    levels_from_root: levels_from_root,
                };
                children.push(level);
            }
        }

        return children;
    }

    /// Adds a node (b) as a sibling of another node (a). They will share a parent node. Distance is
    /// calculated from the parent, not from the sibling.
    pub fn add_sibling(&mut self, a : NodeIndex<u32>, b : String, distance : f32) {
        let parent_node = self.get_parent_node(a);
        self.add_child(parent_node, b, distance);
    }

    /// Adds a child node below another node.
    pub fn add_child(&mut self, parent : NodeIndex<u32>, child : String, distance : f32) {
        let child_node_idx = self.add_node(child);
        self.graph.add_edge(parent, child_node_idx, distance);
    }


    /// Adds a node to the graph.
    pub fn add_node(&mut self, node : String) -> NodeIndex<u32> {
        return self.graph.add_node(node);
    }

    #[allow(dead_code)]
    fn parse_tree_from_string(tree_string: String, graph : &mut Graph<String, f32>, parent : NodeIndex<u32>) -> (&mut Graph<String, f32>, usize) {
        let mut current_leaf = String::new();
        let mut current_distance = String::new();
        let mut branches : Vec<NodeIndex<u32>> = Vec::new();
        let mut branch_distances : Vec<f32> = Vec::new();

        let mut read_mode = "LEAF"; // The read mode can be LEAF, DIST, BDONE

        let mut c = 0;
        let mut chr: char;

        while c < tree_string.len() {
            chr = tree_string.chars().nth(c).expect("Could not read cth character.");
            c += 1;

            if chr == '(' {
                let remainder = tree_string.chars().skip(c).collect();
                let node_count_name = String::from(format!(">>{}", graph.node_count()));

                let branch_node = graph.add_node(node_count_name);
                branches.push(branch_node );
                branch_distances.push(f32::NAN);

                let (_graph, new_skip) = Tree::parse_tree_from_string(
                    remainder, graph, branches[branches.len() - 1]
                );

                c += new_skip;
                read_mode = "BDON";
                continue;
            } else if chr == ')' {
                if read_mode == "DIST" {
                    if current_leaf.len() > 0 {
                        // We found a leaf with a distance.
                        let n = graph.add_node(current_leaf.clone().trim().to_string());
                        graph.add_edge(parent, n, current_distance.parse().unwrap());
                    } else {
                        // We found the distance of a branch, replace NAN with true distance
                        let bdsize = branch_distances.len();
                        branch_distances[bdsize - 1] = current_distance.parse::<f32>().unwrap();
                    }
                } else if read_mode == "LEAF" {
                    // We found a leaf without a distance.
                    let n = graph.add_node(current_leaf.clone().trim().to_string());
                    graph.add_edge(parent, n, f32::NAN);
                }
                break;
            } else if chr == ':' {
                read_mode = "DIST";
            } else if chr == ',' {
                // End of an item in a list
                if read_mode == "DIST" {
                    if current_leaf.len() > 0 {
                        // We found a leaf with a distance.
                        let n = graph.add_node(current_leaf.clone().trim().to_string());
                        graph.add_edge(parent, n, current_distance.parse().unwrap());
                    } else {
                        // We found the distance of a branch, replace NAN with true distance
                        let bdsize = branch_distances.len();
                        branch_distances[bdsize - 1] = current_distance.parse::<f32>().unwrap();
                    }
                } else if read_mode == "LEAF" {
                    // We found a leaf without a distance.
                    let n = graph.add_node(current_leaf.clone().trim().to_string());
                    graph.add_edge(parent, n, f32::NAN);
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

        for (branch, dist) in branches.iter().zip(branch_distances.iter()) {
            graph.add_edge(parent, *branch, *dist);
        }

        return (graph, c)

    }

    /// Takes a newick tree description string and parses it into a Tree instance.
    #[allow(dead_code)]
    pub fn parse(tree_string: String) -> Tree {
        let mut graph : Graph<String, f32> = Graph::<String, f32>::new();
        let parent = graph.add_node(String::from("<root>"));
        Tree::parse_tree_from_string(tree_string, &mut graph, parent);

        let tree = Tree { graph };
        return tree;
    }

    /// Returns a new, empty Tree.
    #[allow(dead_code)]
    pub fn new() -> Tree {
        let graph : Graph<String, f32> = Graph::new();
        Tree { graph }
    }

    /// Finds the index of a node by its name.
    pub fn find_node_idx(&self, name : &String) -> Option<NodeIndex<u32>>{
        for n in self.graph.node_indices() {
            if name == self.graph.node_weight(n).expect("Graph is inconsistent.") {
                return Some(n);
            }
        }
        return None;
    }

    /// Finds the number of levels a node is removed from the root node.
    pub fn get_levels_from_root(&self, n : NodeIndex<u32>) -> usize {
        let mut current_node = n;
        let mut levels_from_root = 0;
        while self.graph.node_weight(current_node).expect("No label for node.") != "<root>" {
            levels_from_root += 1;
            current_node = self.get_parent_node(current_node);
        }

        if levels_from_root == 0 {
            return 0;
        } else {
            return levels_from_root - 1;
        }

    }

    /// Converts this Tree into a distance matrix and returns the resulting TreeDistanceMatrix object.
    pub fn to_distance_matrix(&self) -> TreeDistanceMatrix {
        let mut n_leaves: usize = 0;
        let mut max_depth = 0;
        let mut leaf_map: HashMap<String, usize> = HashMap::new();

        for child in self.traverse_children() {
            for leaf in child.leaves.iter() {
                leaf_map.insert(leaf.clone(), n_leaves);
                n_leaves += 1;
            }
            let levels_from_root = self.get_levels_from_root(child.leaf_nodes[0]);
            if levels_from_root > max_depth {
                max_depth = levels_from_root;
            }
        }

        let mut leaf_distance_matrix: Array2<f32> = Array2::zeros((n_leaves, max_depth + 1));
        let mut identity_matrix: Array2<usize> = Array2::zeros((n_leaves, max_depth + 1));

        let mut finished_leaves: Vec<usize> = Vec::new();
        let mut accumulated_distances: Vec<f32> = Vec::new();
        let mut internal_nodes: Vec<usize> = Vec::new();

        let mut previous_level = 0;
        let mut current_level;
        let mut internal_node_count = 0;

        for (_c, level) in self.traverse_children().iter().enumerate() {
            current_level = level.levels_from_root;

            if previous_level < current_level {
                for _l in current_level..previous_level {
                    accumulated_distances.pop();
                    internal_nodes.pop();
                }
            }
            if previous_level == current_level {
                // Switch branches on the same level.
                accumulated_distances.pop();
                internal_nodes.pop();
            }
            // Compensate for large jumps,
            if (current_level as i32 - previous_level as i32) > 1 {
                for _l in previous_level..(current_level - 1) {
                    accumulated_distances.push(0.0);
                    internal_nodes.push(internal_node_count);
                    internal_node_count += 1;
                }
            }
            if (level.level_distance).is_nan() {
                accumulated_distances.push(0.0);
            } else {
                accumulated_distances.push(level.level_distance);
            }

            internal_nodes.push(internal_node_count);
            internal_node_count += 1;

            for (i, leaf) in level.leaves.iter().enumerate() {
                let leaf_id = leaf_map[leaf];
                for l in 0..current_level {
                    let sum_dist = accumulated_distances[0..l + 1].iter().sum();
                    leaf_distance_matrix[[leaf_id, l]] = sum_dist;
                    identity_matrix[[leaf_id, l]] = internal_nodes[l];
                }
                // Add the final accumulated leaf distance
                leaf_distance_matrix[[leaf_id, current_level]] = leaf_distance_matrix[[leaf_id, current_level - 1]] + level.leaf_distances[i];
                finished_leaves.push(leaf_id);
            }

            previous_level = current_level;
        }

        let distance_matrix = TreeDistanceMatrix::new(
            leaf_distance_matrix, identity_matrix, leaf_map
        );

        distance_matrix
    }

    /// Converts the tree into a dot string that can be visualized with graph-viz.
    pub fn to_dot(&self) -> String {
        let dot_string = Dot::with_config(&self.graph, &[Config::EdgeNoLabel]);
        return format!("{:?}", dot_string);
    }
}
