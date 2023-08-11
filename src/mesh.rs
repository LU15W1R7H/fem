use crate::{fx, geo::Tria};

pub struct TriaMesh2d {
  pub vertex_coords: na::Matrix2xX<fx>,
  pub elements: na::Matrix3xX<usize>,
}

impl TriaMesh2d {
  pub fn tria(&self, global_idx: usize) -> Tria {
    let mut vertex_coords = na::Matrix2x3::zeros();
    for i in 0..3 {
      vertex_coords
        .column_mut(i)
        .copy_from(&self.vertex_coords.column(self.elements[(i, global_idx)]));
    }
    Tria { vertex_coords }
  }
}

pub fn load_tria_mesh_2d(path: impl AsRef<std::path::Path>) -> TriaMesh2d {
  use std::{
    fs,
    io::{self, BufRead},
  };
  let file = fs::File::open(path).unwrap();
  let reader = io::BufReader::new(file);
  let mut lines = reader.lines();

  // read vertices
  let first_line = lines.next().unwrap().unwrap();
  let parts: Vec<&str> = first_line.split_whitespace().collect();
  let nvertices: usize = parts[0].parse().unwrap();
  assert_eq!(parts[1], "Vertices");
  let mut vertex_coords = na::Matrix2xX::zeros(nvertices);
  for i in 0..nvertices {
    let line = lines.next().unwrap().unwrap();
    let parts: Vec<fx> = line
      .split_whitespace()
      .map(|s| s.parse().unwrap())
      .collect();
    vertex_coords[(0, i)] = parts[0];
    vertex_coords[(1, i)] = parts[1];
  }

  // read elements
  let second_header = lines.next().unwrap().unwrap();
  let parts: Vec<&str> = second_header.split_whitespace().collect();
  let nelements: usize = parts[0].parse().unwrap();
  let mut elements = na::Matrix3xX::zeros(nelements);
  for i in 0..nelements {
    let line = lines.next().unwrap().unwrap();
    let parts: Vec<usize> = line
      .split_whitespace()
      .map(|s| s.parse().expect("Expected an integer"))
      .collect();
    elements[(0, i)] = parts[0];
    elements[(1, i)] = parts[1];
    elements[(2, i)] = parts[2];
  }

  TriaMesh2d {
    vertex_coords,
    elements,
  }
}
