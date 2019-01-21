use std::io::BufReader;
use std::io::prelude::*;
use std::fs::File;

use std::collections::BTreeMap;


pub fn batch_id_tax_mapping(identifiers : Vec<String>, tax_id_map : &BTreeMap<[u8; 10], u32>) -> Vec<u32> {
    let u8_identifiers : Vec<[u8; 10]> = identifiers.into_iter().map(uniprot_id_to_u8_arr).collect();
    let result : Vec<u32> = u8_identifiers.iter().map(|id| get_from_map(id, &tax_id_map)).collect();
    return result;
}

pub fn get_from_map(id : &[u8; 10], tax_id_map : &BTreeMap<[u8; 10], u32>) -> u32 {
    match tax_id_map.get(id) {
        Some(i) => { return *i; }
        None => { return 0; }
    }
}

pub fn uniprot_id_to_u8_arr(source_string : String) -> [u8; 10]{
    let mut uniprot_id : [u8; 10] = [0 as u8; 10];
    for (b, byte) in source_string.bytes().enumerate() {
        if b == 10 {
            break;
        }
        uniprot_id[b] = byte;
    }
    return uniprot_id;
}


pub fn load_mapping_file(filename : &String) -> BTreeMap<[u8; 10], u32> {
    let f = File::open(filename).expect("file not found");
    let mut f = BufReader::new(f);

    let mut tax_id_map : BTreeMap<[u8; 10], u32> = BTreeMap::new();

    let mut line = String::new();

    while f.read_line(&mut line).unwrap() > 0 {
        // Wrap this in an extra block to allow the immutable borrow
        {
            let ls : Vec<&str> = line.split_whitespace().collect();

            if ls.len() < 2 {
                continue;
            }

            let uniprot_id = uniprot_id_to_u8_arr(ls[0].to_string());
            let tax_id = ls[2].parse::<u32>().unwrap();

            tax_id_map.insert(uniprot_id, tax_id);
        }
        line.clear();
    }

    return tax_id_map;
}