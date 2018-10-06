extern crate ndarray;

mod tree;
mod utils;
mod tests;
mod test_tree_parsing;


fn main() {
    let tree_file = utils::load_tree_file(String::from(
        "/home/stephan/newick_trees/1.tree"
    ));
    let parsed_tree = tree::Tree::parse(tree_file);
    let children = parsed_tree.traverse_children();
    println!("{:?}", children)
}