use std::collections::BTreeMap;
use std::collections::HashMap;
use std::collections::HashSet;

use std::error::Error;
use std::fs::File;
use std::io::prelude::*;
use std::path::Path;
use std::time;

use bincode::{serialize, deserialize};
use graph_tree;
use tree_distance_matrix::TreeDistanceMatrix;


pub fn load_tree_file(filename: String) -> String {
    let mut f = File::open(filename).expect("file not found");

    let mut contents = String::new();
    f.read_to_string(&mut contents)
        .expect("something went wrong reading the file");

    // Strip newlines from the file.
    String::from(contents.trim())
}

#[allow(dead_code)]
pub fn keys_vec<T: std::cmp::Eq + std::hash::Hash + std::clone::Clone, V>(
    hm: &HashMap<T, V>,
) -> Vec<T> {
    return hm.keys().cloned().collect();
}

#[allow(dead_code)]
pub fn vec_to_set<T: std::cmp::Eq + std::hash::Hash + std::clone::Clone>(
    vector: &Vec<T>,
) -> HashSet<T> {
    let mut set: HashSet<T> = HashSet::new();
    for element in vector.iter() {
        set.insert(element.clone());
    }
    return set;
}

#[allow(dead_code)]
pub fn argsort<T: std::cmp::PartialOrd + std::clone::Clone>(sort_vec : &Vec<T>) -> Vec<usize> {
    let mut enum_vec = sort_vec.iter().enumerate().collect::<Vec<(usize, &T)>>();
    enum_vec.sort_by(|a, &b| a.1.partial_cmp(b.1).unwrap());
    return enum_vec.iter().map(|x| x.0).collect::<Vec<usize>>();
}

pub fn get_indices_vec<T: std::cmp::Eq + std::clone::Clone>(vector : &Vec<T>, value : T ) -> Vec<usize> {
    vector.iter().enumerate().filter(| e | *e.1 == value).map(| e | e.0).collect::<Vec<usize>>()
}

/* Caching utilities */
pub fn cache_mapping(mapping : &BTreeMap<[u8; 10], u32>, path : &String) -> Result<(), Box<Error>> {
    let serialized = serialize(mapping).unwrap();
    let mut file = File::create(&path)?;
    file.write_all(&serialized)?;
    Ok(())
}

pub fn load_cache_mapping(path : &String) -> Result<BTreeMap<[u8; 10], u32>, Box<Error>> {
    let mut file = File::open(path)?;
    let mut contents : Vec<u8> = Vec::new();
    file.read_to_end(&mut contents)?;
    let data: BTreeMap<[u8; 10], u32> = deserialize(&contents)?;
    return Ok(data);
}

pub fn get_cache_mapping_name(original_path : &String) -> String {
    return format!("cache/{}.bin", Path::new(original_path).file_name().unwrap().to_str().unwrap())
}

pub fn convert_file_to_distance_matrix(fname : String) -> TreeDistanceMatrix {
    let now = time::Instant::now();
    let tree_file = load_tree_file(String::from(fname));
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

    let now = time::Instant::now();
    let distance_matrix = parsed_tree.to_distance_matrix();
    println!(
        "Built distance_matrix in {}.{} seconds.",
        now.elapsed().as_secs(),
        now.elapsed().subsec_millis()
    );

    return distance_matrix;
}