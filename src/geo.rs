use crate::fx;

pub struct Tria {
  pub vertex_coords: na::Matrix2x3<fx>,
}

impl Tria {
  pub fn area(&self) -> fx {
    let v = &self.vertex_coords;
    (0.5
      * ((v[(0, 1)] - v[(0, 0)]) * (v[(1, 2)] - v[(1, 1)])
        - (v[(0, 2)] - v[(0, 1)]) * (v[(1, 1)] - v[(1, 0)])))
      .abs()
  }
}
