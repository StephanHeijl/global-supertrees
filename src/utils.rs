use std::fs::File;
use std::io::prelude::*;
use std::collections::HashMap;
use std::collections::HashSet;


pub fn load_tree_file(filename : String) -> String {
    let mut f = File::open(filename).expect("file not found");

    let mut contents = String::new();
    f.read_to_string(&mut contents)
        .expect("something went wrong reading the file");

    // Strip newlines from the file.
    String::from(contents.trim())
}

pub fn keys_vec<T : std::cmp::Eq + std::hash::Hash + std::clone::Clone, V>(hm : &HashMap<T, V>) -> Vec<T> {
    return hm.keys().cloned().collect();
}

pub fn vec_to_set<T : std::cmp::Eq + std::hash::Hash + std::clone::Clone>(vector : &Vec<T>) -> HashSet<T> {
    let mut set : HashSet<T> = HashSet::new();
    for element in vector.iter() {
        set.insert(element.clone());
    }
    return set;
}