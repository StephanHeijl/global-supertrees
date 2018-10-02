#[cfg(test)]
mod tests {
    use::tree;
    use::utils;

    #[test]
    fn test_construct_tree_from_file_1() {
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

    #[test]
    fn test_construct_tree_from_file_2() {
        // Tree 2
        let filename = String::from("/home/stephan/newick_trees/2.tree");
        let tree_file = utils::load_tree_file(filename);
        let parsed_tree = tree::Tree::parse(tree_file);

        let expected_children = vec!(
            vec!(String::from("Bovine"), String::from("Mouse")),
            vec!(String::from("Gibbon")),
            vec!(String::from("Orang")),
            vec!(String::from("Gorilla")),
            vec!(String::from("Chimp"), String::from("Human")),
        );
        let children = parsed_tree.traverse_children();

        // Depth first check
        for (c, child) in children.iter().enumerate() {
            assert_eq!(child[0].0.leaves, expected_children[c]);
        }
    }

    #[test]
    fn test_construct_tree_from_file_3() {
        let filename = String::from("/home/stephan/newick_trees/4.tree");
        let tree_file = utils::load_tree_file(filename);
        let parsed_tree = tree::Tree::parse(tree_file);

        let expected_children = vec!(
            vec!(
                String::from("Alpha"),
                String::from("Beta"),
                String::from("Gamma"),
                String::from("Delta"),
                String::from(""),
                String::from("Epsilon"),
                String::from(""),
                String::from(""),
                String::from("")
            )
        );
        let children = parsed_tree.traverse_children();

        // Depth first check
        for (c, child) in children.iter().enumerate() {
            assert_eq!(child[0].0.leaves, expected_children[c]);
        }
    }
}