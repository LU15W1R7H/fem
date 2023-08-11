use crate::{fx, geo::Tria};

pub fn barycentric_gradients(tria: &Tria) -> na::Matrix2x3<fx> {
  let mut x = na::Matrix3::zeros();
  x.fixed_view_mut::<3, 1>(0, 0)
    .copy_from(&na::Vector3::from_element(1.0));
  x.fixed_view_mut::<3, 2>(0, 1)
    .copy_from(&tria.vertex_coords.transpose());
  x.try_inverse().unwrap().fixed_view::<2, 3>(1, 0).into()
}

#[rustfmt::skip]
pub fn elmat_mass(tria: &Tria) -> na::Matrix3<fx> {
  na::Matrix3::new(
    2.0, 1.0, 1.0,
    1.0, 2.0, 1.0,
    1.0, 1.0, 2.0
  ) * tria.area()
    / 12.0
}

pub fn elmat_laplacian(tria: &Tria) -> na::Matrix3<fx> {
  let x = barycentric_gradients(tria);
  tria.area() * x.transpose() * x
}

pub fn elvec<F>(tria: &Tria, f: F) -> na::Vector3<fx>
where
  F: Fn(na::VectorView2<fx>) -> fx,
{
  let mut phi = na::Vector3::zeros();
  for i in 0..3 {
    phi[i] = f(tria.vertex_coords.column(i));
  }
  phi *= tria.area() / 3.0;
  phi
}
