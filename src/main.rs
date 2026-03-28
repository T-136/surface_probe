use std::collections::HashMap;

use surface_probe::{
    occ_onlyocc_from_xyz, read_and_write::AtomPosition, read_atom_sites, read_nn, read_sample,
    surface, write_occ_as_xyz,
};

fn main() {
    let inp = "../../kmc_sims/201A_coat_900K_1.xyz";
    let grid_folder = "../kmc_cluster_simulation_rs/303030-grid_kmc";

    let input_cluster = read_sample(inp);
    let nn = read_nn(&format!("{}/nearest_neighbor", grid_folder));
    let xsites_positions = read_atom_sites(&format!("{}/atom_sites", grid_folder), 0);
    let mut atom_pos: Vec<AtomPosition> = vec![
        AtomPosition {
            occ: 255,
            ..Default::default()
        };
        xsites_positions.len()
    ];
    let mut atom_names: HashMap<String, u8> = HashMap::new();
    atom_names.insert("Pt".to_string(), 0);
    atom_names.insert("Pd".to_string(), 1);
    let onlyocc = occ_onlyocc_from_xyz(
        &mut atom_pos,
        &input_cluster,
        &xsites_positions,
        &atom_names,
        None,
        &nn,
    );
    let (surface, empty) = surface(&onlyocc, &nn, &xsites_positions);
    let mut occ_surface = vec![255; atom_pos.len()];
    occ_surface
        .iter_mut()
        .enumerate()
        .filter(|x| surface.contains(&(x.0 as u32)))
        .for_each(|(i, x)| *x = atom_pos[i].occ);
    let mut occ_near_surface = vec![255; atom_pos.len()];
    occ_near_surface
        .iter_mut()
        .enumerate()
        .filter(|x| empty.contains(&(x.0 as u32)))
        .for_each(|(i, x)| *x = 1);

    let v = vec![occ_surface, occ_near_surface];

    write_occ_as_xyz(
        "surface.xyz",
        "test_output".to_string(),
        &v,
        &xsites_positions,
        &[0.; 3],
        &atom_names,
    );
}
