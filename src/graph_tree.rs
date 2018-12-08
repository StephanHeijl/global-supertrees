use ndarray::prelude::*;
use std::collections::HashMap;
use std::f32;
use tree_distance_matrix::*;
use petgraph::prelude::NodeIndex;
use petgraph::visit::Dfs;
use petgraph::{Graph, Incoming};

#[derive(Debug)]
pub struct Tree {
    pub graph : Graph<String, f32>,
}


impl Tree {

    pub fn traverse_children(&self) -> Vec<Vec<String>> {

        let mut children : Vec<Vec<String>> = Vec::new();
        let mut current_branch : Vec<String> = Vec::new();
        let mut dfs = Dfs::new(&self.graph, self.graph.node_indices().nth(0).expect(""));

        while let Some(node) = dfs.next(&self.graph) {
            let node_name = self.graph.node_weight(node).unwrap();
            if node_name.starts_with(">>") {
                children.push(current_branch.clone());
                current_branch.clear();
            } else {
                current_branch.push(node_name.to_string());
            }
        }

        children.remove(0);
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

                let (graph, new_skip) = Tree::parse_tree_from_string(
                    remainder,  graph, branches[branches.len() - 1]
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

        Tree{
            graph
        }
    }

    #[allow(dead_code)]
    pub fn new() -> Tree {
        let graph : Graph<String, f32> = Graph::new();

        Tree { graph }
    }
}
