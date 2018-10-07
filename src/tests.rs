#[cfg(test)]
mod tests {
    use::tree;
    use::utils;

    #[test]
    fn test_traversal() {
        let mut tree = tree::Tree::new (
            vec!(String::from("Bovine"), String::from("Wolf")),
            vec!(),
            vec!(0.1, 0.4),
            vec!()
        );

        tree.add_leaf(String::from("Homo sapiens"), 0.6);

        let mut tree_two = tree::Tree::new (
            vec!(String::from("T2")),
            vec!(),
            vec!(0.1),
            vec!()
        );
        let tree_three = tree::Tree::new (
            vec!(String::from("T3")),
            vec!(),
            vec!(0.2),
            vec!()
        );
        let tree_four = tree::Tree::new (
            vec!(String::from("T4")),
            vec!(),
            vec!(0.3),
            vec!()
        );

        tree_two.add_branch(tree_four, 2.0);
        tree.add_branch(tree_two, 10.0);
        tree.add_branch(tree_three, 2.0);

        let tc = tree.traverse_children();
        let expected_leaves = vec!("T2", "T4", "T3");
        for (i, hm) in tc.iter().enumerate() {
            assert_eq!(hm[0].0.leaves[0], expected_leaves[i]);
        }
    }

    #[test]
    fn test_load_tree_file() {
        let filename = String::from("/home/stephan/newick_trees/1.tree");
        let tree_file = utils::load_tree_file(filename);
        //println!("{}" ,tree_file);
        assert!(tree_file.contains("("));
    }

    #[test]
    fn test_to_distance_matrix() {
        let filename = String::from("/home/stephan/newick_trees/1.tree");
        let tree_file = utils::load_tree_file(filename);
        let parsed_tree = tree::Tree::parse(tree_file);
        //println!("{:?}", parsed_tree);
        parsed_tree.to_distance_matrix();
    }
}