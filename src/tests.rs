#[cfg(test)]
mod tests {
    use ndarray::prelude::*;
    use tree;
    use tree_distance_matrix;
    use tree_merging;
    use utils;

    #[test]
    fn test_traversal() {
        let mut tree = tree::Tree::new(
            vec![String::from("Bovine"), String::from("Wolf")],
            vec![],
            vec![0.1, 0.4],
            vec![],
        );

        tree.add_leaf(String::from("Homo sapiens"), 0.6);

        let mut tree_two = tree::Tree::new(vec![String::from("T2")], vec![], vec![0.1], vec![]);
        let tree_three = tree::Tree::new(vec![String::from("T3")], vec![], vec![0.2], vec![]);
        let tree_four = tree::Tree::new(vec![String::from("T4")], vec![], vec![0.3], vec![]);

        tree_two.add_branch(tree_four, 2.0);
        tree.add_branch(tree_two, 10.0);
        tree.add_branch(tree_three, 2.0);

        let tc = tree.traverse_children();
        let expected_leaves = vec!["T2", "T4", "T3"];
        for (i, hm) in tc.iter().enumerate() {
            assert_eq!(hm.0.leaves[0], expected_leaves[i]);
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
    fn test_find_first_zero() {
        let arr_one = Array::from_vec(vec![0, 1, 2, 3, 4, 5]);
        let arr_two = Array::from_vec(vec![0, 1, 2, 0, 4, 5]);
        let arr_three = Array::from_vec(vec![0, 1, 2, 3, 4, 0]);
        let arr_four = Array::from_vec(vec![1, 2, 3, 4, 5, 6]);
        assert_eq!(
            tree_distance_matrix::TreeDistanceMatrix::find_first_zero(arr_one.view()),
            0
        );
        assert_eq!(
            tree_distance_matrix::TreeDistanceMatrix::find_first_zero(arr_two.view()),
            3
        );
        assert_eq!(
            tree_distance_matrix::TreeDistanceMatrix::find_first_zero(arr_three.view()),
            5
        );
        assert_eq!(
            tree_distance_matrix::TreeDistanceMatrix::find_first_zero(arr_four.view()),
            0
        );
    }

    #[test]
    fn test_find_first_common_ancestor() {
        let identity_matrix = Array::from_shape_vec(
            (4, 5),
            vec![0, 1, 2, 3, 4, 0, 1, 2, 5, 6, 0, 7, 8, 9, 0, 0, 7, 8, 9, 10],
        )
        .unwrap();
        assert_eq!(
            tree_distance_matrix::TreeDistanceMatrix::find_first_common_ancestor(
                0,
                1,
                &identity_matrix
            ),
            2
        );
        assert_eq!(
            tree_distance_matrix::TreeDistanceMatrix::find_first_common_ancestor(
                1,
                0,
                &identity_matrix
            ),
            2
        );
        assert_eq!(
            tree_distance_matrix::TreeDistanceMatrix::find_first_common_ancestor(
                1,
                2,
                &identity_matrix
            ),
            0
        );
        assert_eq!(
            tree_distance_matrix::TreeDistanceMatrix::find_first_common_ancestor(
                2,
                1,
                &identity_matrix
            ),
            0
        );
        assert_eq!(
            tree_distance_matrix::TreeDistanceMatrix::find_first_common_ancestor(
                2,
                3,
                &identity_matrix
            ),
            3
        );
        assert_eq!(
            tree_distance_matrix::TreeDistanceMatrix::find_first_common_ancestor(
                3,
                2,
                &identity_matrix
            ),
            3
        );
        assert_eq!(
            tree_distance_matrix::TreeDistanceMatrix::find_first_common_ancestor(
                0,
                2,
                &identity_matrix
            ),
            0
        );
        assert_eq!(
            tree_distance_matrix::TreeDistanceMatrix::find_first_common_ancestor(
                2,
                0,
                &identity_matrix
            ),
            0
        );
    }

    #[test]
    fn test_to_distance_matrix_simple() {
        let filename = String::from("/home/stephan/newick_trees/5.tree");
        let tree_file = utils::load_tree_file(filename);
        let parsed_tree = tree::Tree::parse(tree_file);
        parsed_tree.to_distance_matrix();
    }

    #[test]
    fn test_to_distance_matrix_hard() {
        let filename = String::from("/home/stephan/newick_trees/1.tree");
        let tree_file = utils::load_tree_file(filename);
        let parsed_tree = tree::Tree::parse(tree_file);
        parsed_tree.to_distance_matrix();
    }

    #[test]
    fn test_neighbour_joining() {
        let tree_string = String::from("(Bovine:0.69395,(Hylobates:0.36079,(Pongo:0.33636,(G._Gorilla:0.17147,(P._paniscus:0.19268,H._sapiens:0.11927):0.08386):0.06124):0.15057):0.54939,Rodent:1.21460);");
        let parsed_tree = tree::Tree::parse(tree_string);
        let distance_matrix = parsed_tree.to_distance_matrix();
        let nj_tree = distance_matrix.neighbour_joining();
        assert_eq!(parsed_tree, nj_tree);
    }

    #[test]
    fn test_tree_merging() {
        let tree_string_one = String::from("(Bovine:0.69395,Lamb:0.2,(Hylobates:0.36079,(Pongo:0.33636,(G._Gorilla:0.17147,(P._paniscus:0.19268,H._sapiens:0.11927):0.08386):0.06124):0.15057):0.54939,Rodent:1.21460);");
        let tree_string_two = String::from("(Bovine:0.29395,Porcine:0.4,(Hylobates:0.16079,(Pongo:0.53636,(G._Gorilla:0.07147,(P._paniscus:0.09268,H._sapiens:0.41927):0.01386):0.12124):0.45057):0.34939,Rodent:2.21460);");

        let tree_one = tree::Tree::parse(tree_string_one);
        let tree_two = tree::Tree::parse(tree_string_two);

        tree_merging::merge_trees(vec![tree_one, tree_two]);
    }

    #[test]
    fn test_get_partial_tree() {
        let tree_string = String::from("(Bovine:0.69395,(Hylobates:0.36079,(Pongo:0.33636,(G._Gorilla:0.17147,(P._paniscus:0.19268,H._sapiens:0.11927):0.08386):0.06124):0.15057):0.54939,Rodent:1.21460);");

        let tree = tree::Tree::parse(tree_string);
        let dm = tree.to_distance_matrix();
        dm.get_partial_distance_matrix(&vec![
            String::from("Hylobates"),
            String::from("Pongo"),
            String::from("Bovine"),
        ]);
    }
}
