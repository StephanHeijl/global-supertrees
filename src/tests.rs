#[cfg(test)]
mod tests {
    use::tree;
    use::utils;

    #[test]
    fn test_traversal() {
        let mut tree = tree::Tree {
            leaves: vec!(String::from("Bovine"), String::from("Wolf")),
            branches: vec!(),
            leaf_distances: vec!(0.1, 0.4),
            branch_distances: vec!()
        };

        tree.add_leaf(String::from("Homo sapiens"), 0.6);

        let mut tree_two = tree::Tree {
            leaves: vec!(String::from("T2")),
            branches: vec!(),
            leaf_distances: vec!(0.1),
            branch_distances: vec!()
        };
        let tree_three = tree::Tree {
            leaves: vec!(String::from("T3")),
            branches: vec!(),
            leaf_distances: vec!(0.2),
            branch_distances: vec!()
        };
        let tree_four = tree::Tree {
            leaves: vec!(String::from("T4")),
            branches: vec!(),
            leaf_distances: vec!(0.3),
            branch_distances: vec!()
        };

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
    fn test_construct_tree_from_file() {
        let filename = String::from("/home/stephan/newick_trees/1.tree");
        let tree_file = utils::load_tree_file(filename);
        let parsed_tree = tree::Tree::parse(tree_file);

        let expected_children = vec!(
            vec!(String::from("dog")),
            vec!(String::from("raccoon"), String::from("bear")),
            vec!(),
            vec!(String::from("sea_lion"), String::from("seal")),
            vec!(String::from("weasel")),
            vec!(String::from("monkey"), String::from("cat"))
        );
        let children = parsed_tree.traverse_children();

        // Depth first check
        for (c, child) in children.iter().enumerate() {
            assert_eq!(child[0].0.leaves, expected_children[c]);
        }
    }
}