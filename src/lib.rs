use std::{
    collections::{HashMap, HashSet},
    mem,
};
pub mod read_and_write;
use fnv::{FnvBuildHasher, FnvHashMap, FnvHashSet};
pub use read_and_write::{
    occ_onlyocc_from_xyz, read_atom_sites, read_nn, read_sample, write_occ_as_xyz,
};

pub fn surface(
    onlyocc: &HashSet<u32, FnvBuildHasher>,
    nn: &HashMap<u32, [u32; 12], FnvBuildHasher>,
    sites: &[[f64; 3]],
) -> (HashSet<u32, FnvBuildHasher>, HashSet<u32, FnvBuildHasher>) {
    // let mut visi: HashMap<u32, [u32; super::CN], FnvBuildHasher> =
    //     FnvHashMap::with_capacity_and_hasher(5400, Default::default());
    let mut visited: HashSet<u32, FnvBuildHasher> = FnvHashSet::with_hasher(Default::default());
    let mut current_positions: HashSet<u32, FnvBuildHasher> =
        FnvHashSet::with_hasher(Default::default());
    let mut surface: HashSet<u32, FnvBuildHasher> = FnvHashSet::with_hasher(Default::default());
    let mut near_surface: HashSet<u32, FnvBuildHasher> =
        FnvHashSet::with_hasher(Default::default());
    get_start_positions(sites, &mut current_positions);
    println!("onlyocc: {:?}", onlyocc);

    let mut future_positions: HashSet<u32, FnvBuildHasher> =
        FnvHashSet::with_hasher(Default::default());

    while !current_positions.is_empty() {
        println!("len visited: {}", visited.len());
        for site in current_positions.iter() {
            let neighbors = nn.get(&site).unwrap();
            for neighbor in neighbors {
                if visited.contains(neighbor) || current_positions.contains(neighbor) {
                    continue;
                } else if onlyocc.contains(neighbor) {
                    surface.insert(*neighbor);
                    near_surface.insert(*site);
                } else {
                    future_positions.insert(*neighbor);
                }
            }
        }
        visited.extend(current_positions.iter());
        mem::swap(&mut future_positions, &mut current_positions);

        future_positions.clear();
    }
    (surface, near_surface)
}

pub fn get_start_positions(xyz: &[[f64; 3]], start_sites: &mut HashSet<u32, FnvBuildHasher>) {
    let min_max: [(f64, f64); 3] = [
        min_max_xyz(xyz, 0),
        min_max_xyz(xyz, 1),
        min_max_xyz(xyz, 2),
    ];

    for (grid_index, pos) in xyz.iter().enumerate() {
        for (i, x) in pos.iter().enumerate() {
            if min_max[i].0 == *x || min_max[i].1 == *x {
                start_sites.insert(grid_index as u32);
                break;
            }
        }
    }
}

//reused from kmc
fn min_max_xyz(xyz: &[[f64; 3]], i: usize) -> (f64, f64) {
    let min = xyz
        .iter()
        .map(|p| p[i])
        .min_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap();
    let max = xyz
        .iter()
        .map(|p| p[i])
        .max_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap();
    (min, max)
}
