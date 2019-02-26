#[cfg(test)]
mod tests {

    use ndarray::prelude::*;
    use rayon::prelude::*;
    use graph_tree as tree;
    use tree_distance_matrix;
    use tree_distance_matrix::TreeDistanceMatrix;
    use tree_merging;
    use utils;
    use uniprot;
    use std::collections::BTreeMap;

    #[test]
    fn test_load_tree_file() {
        let filename = String::from("newick_trees/1.tree");
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
    fn test_get_levels_from_root() {
        let filename = String::from("newick_trees/1.tree");
        let tree_file = utils::load_tree_file(filename);
        let parsed_tree = tree::Tree::parse(tree_file);

        let expected_levels : Vec<usize> = vec!(1, 4, 4, 2);
        let nodes : Vec<String> = vec!(
            String::from("dog"),
            String::from("monkey"),
            String::from("cat"),
            String::from("raccoon")
        );

        for (node_name, el) in nodes.iter().zip(expected_levels.iter()) {
            let target = parsed_tree.find_node_idx(node_name).unwrap();
            assert_eq!(el, &parsed_tree.get_levels_from_root(target))
        }
    }

    #[test]
    fn test_to_distance_matrix_simple() {
        let filename = String::from("newick_trees/5.tree");
        let tree_file = utils::load_tree_file(filename);
        let parsed_tree = tree::Tree::parse(tree_file);
        let dm = parsed_tree.to_distance_matrix();
        assert_eq!(dm.get_distance(String::from("Alpha"), String::from("Gamma")), Some(18.));
        assert_eq!(dm.get_distance(String::from("Alpha"), String::from("Delta")), Some(16.));
        assert_eq!(dm.get_distance(String::from("Delta"), String::from("Beta")), Some(30.));
        assert_eq!(dm.get_distance(String::from("Gamma"), String::from("Beta")), Some(32.));
        assert_eq!(dm.get_distance(String::from("Alpha"), String::from("Alpha")), Some(0.));
        println!("{}", dm.to_csv());
    }

    #[test]
    fn test_to_distance_matrix_medium() {
        let filename = String::from("newick_trees/1.tree");
        let tree_file = utils::load_tree_file(filename);
        let parsed_tree = tree::Tree::parse(tree_file);
        let dm = parsed_tree.to_distance_matrix();

        assert_eq!(dm.get_distance(String::from("raccoon"), String::from("bear")), Some(26.000002));
        assert_eq!(dm.get_distance(String::from("raccoon"), String::from("dog")), Some(45.507133));
        println!("{}", dm.to_csv());
    }

    #[test]
    fn test_to_distance_matrix_medium_nested() {
        let filename = String::from("newick_trees/6.tree");
        let tree_file = utils::load_tree_file(filename);
        let parsed_tree = tree::Tree::parse(tree_file);
        let dm = parsed_tree.to_distance_matrix();

        assert_eq!(dm.get_distance(String::from("raccoon"), String::from("bear")), Some(26.000002));
        assert_eq!(dm.get_distance(String::from("raccoon"), String::from("dog")), Some(45.507133));
        println!("{}", dm.to_csv());
    }

    #[test]
    fn test_neighbour_joining() {
        let tree_string = String::from("(Bovine:0.69395,(Hylobates:0.36079,(Pongo:0.33636,(G._Gorilla:0.17147,(P._paniscus:0.19268,H._sapiens:0.11927):0.08386):0.06124):0.15057):0.54939,Rodent:1.21460);");
        let parsed_tree = tree::Tree::parse(tree_string);
        let distance_matrix = parsed_tree.to_distance_matrix();
        println!("{}", distance_matrix.to_csv());

        let nj_tree = distance_matrix.neighbour_joining();
        assert_eq!(parsed_tree, nj_tree);

        println!("{:?}", distance_matrix);
        println!("{:?}", nj_tree.to_distance_matrix());
    }


    #[test]
    fn test_tree_merging() {
        let tree_string_one = String::from("(Bovine:0.69395,Lamb:0.2,Porcine:0.4,Memine:0.5,Zoobelflurb:0.1,(Hylobates:0.36079,(Pongo:0.33636,(G._Gorilla:0.17147,(P._paniscus:0.19268,H._sapiens:0.11927):0.08386):0.06124):0.15057):0.54939,Rodent:1.21460);");
        let tree_string_two = String::from("(Bovine:0.29395,Deer:0.1,(Hylobates:0.16079,Balohytes:0.15,(Pongo:0.53636,(G._Gorilla:0.07147,(P._paniscus:0.09268,H._sapiens:0.41927):0.01386):0.12124):0.45057):0.34939,Rodent:2.21460);");

        let tree_one = tree::Tree::parse(tree_string_one);
        let tree_two = tree::Tree::parse(tree_string_two);

        let merged_tree = tree_merging::merge_trees(
            vec!(tree_one, tree_two)
        );

        println!("{}", merged_tree.to_dot());


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

    #[test]
    fn test_to_dot() {
        let tree_string = String::from("(Bovine:0.69395,(Hylobates:0.36079,(Pongo:0.33636,(G._Gorilla:0.17147,(P._paniscus:0.19268,H._sapiens:0.11927):0.08386):0.06124):0.15057):0.54939,Rodent:1.21460);");
        let tree = tree::Tree::parse(tree_string);
        println!("{}", tree.to_dot());
    }

    #[test]
    fn test_batch_id_tax_mapping() {
        let identifiers = vec!("Q6GZX4", "Q91G88", "Q6GZX2");
        let identifiers = identifiers.iter().map(|s| String::from(*s)).collect();
        let tax_id_map = uniprot::load_mapping_file(&String::from("taxids/taxids_small.txt"));
        println!("{:?}", uniprot::batch_id_tax_mapping(identifiers, &tax_id_map));
    }

    #[test]
    fn test_merge_organisms() {
        let tax_id_map = uniprot::load_mapping_file(&String::from("taxids/taxids_small.txt"));
        let leaves =  vec!("Q197F5", "Q197F7", "Q197F8", "Q6GZW9", "Q6GZX0", "Q6GZX1", "Q6GZX2", "Q6GZX3", "Q6GZX4", "Q91G88");
        let leaves: Vec<String> = leaves.iter().map(|x| x.to_string()).collect();
        let distance_matrix = Array2::ones((leaves.len(), leaves.len()));

        let dm = tree_distance_matrix::TreeDistanceMatrix::new_from_matrix_and_leaves(
            distance_matrix,
            leaves
        );

        let new_dm = dm.merge_organisms(&tax_id_map);
        assert!(new_dm.shape()[0] == 3);
        assert!(new_dm.shape()[1] == 3);
    }

    #[test]
    fn test_merge_trees() {
        let mut mapping : BTreeMap<[u8; 10], u32> = BTreeMap::new();
        let cache_mapping_path = utils::get_cache_mapping_name(&"taxids/taxids_small.txt".to_string());
        match utils::load_cache_mapping(&cache_mapping_path) {
            Ok(m) => {
                mapping = m;
                println!("Used the cached the mapping file here: {:?} ", cache_mapping_path);
            },
            Err(_e) => {
                mapping = uniprot::load_mapping_file(&"taxids/taxids_small.txt".to_string());
                utils::cache_mapping(&mapping, &cache_mapping_path).expect(
                    "Could not cache mapping file."
                );
                println!("Cached the mapping file here: {:?} ", cache_mapping_path);
            }
        }

        println!("Loaded tax id map");
        let trees = vec!(
            "newick_trees/benchmark_tree_10k.random.tree",
            "newick_trees/benchmark_tree_20k.random.tree",
            "newick_trees/benchmark_tree_40k.random.tree",
            //"newick_trees/benchmark_tree_80k.random.tree"
        );

        let distance_matrices : Vec<TreeDistanceMatrix> = trees.par_iter().map(
            |t| utils::convert_file_to_distance_matrix(t.to_string())
        ).collect();

        println!("Loaded trees");

        let mut normalized_trees : Vec<tree::Tree> = distance_matrices.par_iter().map(
            |m| m.merge_organisms(&mapping).neighbour_joining()
        ).collect();

        println!("Merged organism distance matrices");

        tree_merging::merge_trees(normalized_trees);

    }

}
