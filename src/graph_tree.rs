use ndarray::prelude::*;
use std::collections::HashMap;
use std::f32;
use tree_distance_matrix::*;
use petgraph::prelude::NodeIndex;
use petgraph::visit::Dfs;
use petgraph::graph::node_index;
use petgraph::{Graph, Incoming};
use petgraph::algo::dijkstra;
use std::process;

#[derive(Debug)]
pub struct Tree {
    pub graph : Graph<String, f32>,
}

#[derive(Debug)]
pub struct Level {
    pub leaves : Vec<String>,
    pub leaf_nodes : Vec<NodeIndex<u32>>,
    pub level_distance : f32,
    pub leaf_distances : Vec<f32>,
    pub levels_from_root: usize,
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

        println!("{:?}", leaves_self);
        println!("{:?}", leaves_other);

        return leaves_self == leaves_other;
    }
}


impl Tree {
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

            match self.graph.edges_directed(node, Incoming).nth(0) {
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
            chr = tree_string.chars().nth(c).unwrap();
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

    #[allow(dead_code)]
    pub fn parse(tree_string: String) -> Tree {
        let mut graph : Graph<String, f32> = Graph::<String, f32>::new();
        let parent = graph.add_node(String::from("<root>"));
        Tree::parse_tree_from_string(tree_string, &mut graph, parent);

        let tree = Tree { graph };

        println!("{:?}", tree.get_leaves());
        return tree;
    }

    #[allow(dead_code)]
    pub fn new() -> Tree {
        let graph : Graph<String, f32> = Graph::new();
        Tree { graph }
    }

    pub fn find_node_idx(&self, name : &String) -> Option<NodeIndex<u32>>{
        for n in self.graph.node_indices() {
            if name == self.graph.node_weight(n).expect("Graph is inconsistent.") {
                return Some(n);
            }
        }
        return None;
    }

    pub fn get_levels_from_root(&self, n : NodeIndex<u32>) -> usize {
        let shortest_path = dijkstra(
            &self.graph,
            node_index(0),
            Some(n),
            | _edge | 1
        );
        let sp = shortest_path[&n];
        //println!("{:?}, {:?} => {}", n, self.graph.node_weight(n), sp);
        if sp == 0 {
            return 0;
        } else {
            return sp - 1;
        }
    }

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

        println!("{:?}", leaf_distance_matrix.shape());
        println!("{:?}", identity_matrix.shape());

        let mut finished_leaves: Vec<usize> = Vec::new();
        let mut accumulated_distances: Vec<f32> = Vec::new();
        let mut internal_nodes: Vec<usize> = Vec::new();

        let mut previous_level = 0;
        let mut current_level;
        let mut internal_node_count = 0;

        for (c, level) in self.traverse_children().iter().enumerate() {
            current_level = level.levels_from_root;

            println!("{:?} ... {} -> {}", level, previous_level, current_level);

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
                for l in previous_level..(current_level - 1) {
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
}
