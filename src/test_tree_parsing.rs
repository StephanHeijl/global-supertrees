#[cfg(test)]
mod tests {
    use tree;

    #[test]
    fn test_construct_tree_from_file_1() {
        let tree_file = String::from("((raccoon:19.19959,bear:6.80041):0.84600,((sea_lion:11.99700, seal:12.00300):7.52973,((monkey:100.85930,cat:47.14069):20.59201, weasel:18.87953):2.09460):3.87382,dog:25.46154);");
        let parsed_tree = tree::Tree::parse(tree_file);

        let expected_children = vec![
            vec![String::from("dog")],
            vec![String::from("raccoon"), String::from("bear")],
            vec![],
            vec![String::from("sea_lion"), String::from("seal")],
            vec![String::from("weasel")],
            vec![String::from("monkey"), String::from("cat")],
        ];
        let children = parsed_tree.traverse_children();

        // Depth first check
        for (c, child) in children.iter().enumerate() {
            assert_eq!(child.0.leaves, expected_children[c]);
        }
    }

    #[test]
    fn test_construct_tree_from_file_2() {
        // Tree 2
        let tree_file = String::from("(Bovine:0.69395,(Gibbon:0.36079,(Orang:0.33636,(Gorilla:0.17147,(Chimp:0.19268, Human:0.11927):0.08386):0.06124):0.15057):0.54939,Mouse:1.21460):0.10;");
        let parsed_tree = tree::Tree::parse(tree_file);

        let expected_children = vec![
            vec![String::from("Bovine"), String::from("Mouse")],
            vec![String::from("Gibbon")],
            vec![String::from("Orang")],
            vec![String::from("Gorilla")],
            vec![String::from("Chimp"), String::from("Human")],
        ];
        let children = parsed_tree.traverse_children();

        // Depth first check
        for (c, child) in children.iter().enumerate() {
            assert_eq!(child.0.leaves, expected_children[c]);
        }
    }

    #[test]
    fn test_construct_tree_from_file_3() {
        let tree_file = String::from("(Alpha,Beta,Gamma,Delta,,Epsilon,,,);");
        let parsed_tree = tree::Tree::parse(tree_file);

        let expected_leaves = vec![
            String::from("Alpha"),
            String::from("Beta"),
            String::from("Gamma"),
            String::from("Delta"),
            String::from(""),
            String::from("Epsilon"),
            String::from(""),
            String::from(""),
            String::from(""),
        ];
        let children = parsed_tree.traverse_children();

        for (i, leaf) in children[0].0.leaves.iter().enumerate() {
            assert_eq!(leaf, &expected_leaves[i]);
        }
    }
}
