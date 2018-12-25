#[macro_use(s)]
extern crate ndarray;
extern crate rayon;
extern crate petgraph;
extern crate regex;

use std::env;
use std::time;

mod test_tree_parsing;
mod tests;
mod graph_tree;
mod tree_distance_matrix;
mod tree_merging;
mod utils;


fn main() {
    let mut fname = String::from("/home/stephan/newick_trees/1.tree");
    for (a, arg) in env::args().enumerate() {
        if a == 1 {
            fname = arg;
        }
    }

    let now = time::Instant::now();
    let tree_file = utils::load_tree_file(String::from(fname));
    println!(
        "Loaded tree in {}.{} seconds.",
        now.elapsed().as_secs(),
        now.elapsed().subsec_millis()
    );

    let now = time::Instant::now();
    let parsed_tree = graph_tree::Tree::parse(tree_file);
    println!(
        "Parsed tree in {}.{} seconds.",
        now.elapsed().as_secs(),
        now.elapsed().subsec_millis()
    );

    let now = time::Instant::now();
    let _children = parsed_tree.traverse_children();
    println!(
        "Built tree traversal map in {}.{} seconds.",
        now.elapsed().as_secs(),
        now.elapsed().subsec_millis()
    );
    //println!("{:?}", children)

    let now = time::Instant::now();
    let distance_matrix = parsed_tree.to_distance_matrix();
    println!(
        "Built distance_matrix in {}.{} seconds.",
        now.elapsed().as_secs(),
        now.elapsed().subsec_millis()
    );
    println!("{:?}", distance_matrix.shape());

    //println!("{}", distance_matrix.to_csv());
}
