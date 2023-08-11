pub use crate::fx;
use crate::{geo::Tria, mesh::TriaMesh2d};

pub fn assemble_galmat<ElmatProvider>(
  mesh: &TriaMesh2d,
  elmat_provider: ElmatProvider,
) -> nas::CscMatrix<fx>
where
  ElmatProvider: Fn(&Tria) -> na::Matrix3<fx>,
{
  let n = mesh.vertex_coords.ncols();
  let m = mesh.elements.ncols();

  let mut galmat = nas::CooMatrix::new(n, n);

  for i in 0..m {
    let tria = mesh.tria(i);
    let elmat = elmat_provider(&tria);
    let element = mesh.elements.column(i);
    for j in 0..3 {
      for k in 0..3 {
        galmat.push(element[j], element[k], elmat[(j, k)]);
      }
    }
  }
  nas::CscMatrix::from(&galmat)
}

pub fn assemble_galvec<ElvecProvider>(
  mesh: &TriaMesh2d,
  elvec_provider: ElvecProvider,
) -> na::DVector<fx>
where
  ElvecProvider: Fn(&Tria) -> na::Vector3<fx>,
{
  let m = mesh.elements.ncols();
  let n = mesh.vertex_coords.ncols();

  let mut galvec = na::DVector::zeros(n);

  for i in 0..m {
    let tria = mesh.tria(i);
    let elvec = elvec_provider(&tria);
    let element = mesh.elements.column(i);
    for j in 0..3 {
      galvec[element[j]] += elvec[j];
    }
  }
  galvec
}
