extern crate nalgebra as na;
extern crate nalgebra_sparse as nas;

use fem::{
  assemble::{assemble_galmat, assemble_galvec},
  fx,
  mesh::{load_tria_mesh_2d, TriaMesh2d},
  uscalfe::{barycentric_gradients, elmat_laplacian, elmat_mass, elvec},
};

fn error_l2<F>(mesh: &TriaMesh2d, mu_fem: na::DVectorView<fx>, u_analytic: F) -> fx
where
  F: Fn(na::VectorView2<fx>) -> fx,
{
  let mut error_sq = 0.0;

  let m = mesh.elements.ncols();
  for i in 0..m {
    let tria = mesh.tria(i);

    let mut error_vertices = na::Vector3::zeros();
    for k in 0..3 {
      error_vertices[k] = u_analytic(tria.vertex_coords.column(k)) - mu_fem[mesh.elements[(k, i)]];
    }
    error_sq += tria.area() / 3.0 * error_vertices.norm_squared();
  }
  error_sq.sqrt()
}

fn error_h1s<F>(mesh: &TriaMesh2d, mu_fem: na::DVectorView<fx>, u_grad_analytic: F) -> fx
where
  F: Fn(na::VectorView2<fx>) -> na::Vector2<fx>,
{
  let mut error_sq = 0.0;

  let m = mesh.elements.ncols();
  for i in 0..m {
    let tria = mesh.tria(i);

    let mut fem_values_vertices = na::Vector3::zeros();
    for k in 0..3 {
      fem_values_vertices[k] = mu_fem[mesh.elements[(k, i)]];
    }
    let fem_gradient_tria = barycentric_gradients(&tria) * fem_values_vertices;

    let mut analytic_gradients_vertices = na::Matrix2x3::zeros();
    for k in 0..3 {
      analytic_gradients_vertices
        .column_mut(k)
        .copy_from(&u_grad_analytic(tria.vertex_coords.column(k)));
    }

    let mut error_vertices = na::Vector3::zeros();
    for k in 0..3 {
      error_vertices[k] =
        (fem_gradient_tria - analytic_gradients_vertices.column(k)).norm_squared();
    }
    error_sq += tria.area() / 3.0 * error_vertices.sum();
  }
  error_sq.sqrt()
}

fn main() {
  use std::f64::consts::{PI, TAU};

  let mesh = load_tria_mesh_2d("res/mesh.txt");
  dbg!(mesh.vertex_coords.ncols());
  dbg!(mesh.elements.ncols());

  let f = |x: na::VectorView2<fx>| -> fx {
    (1.0 + 8.0 * PI * PI) * (TAU * x[0]).cos() * (TAU * x[1]).cos()
  };

  let u_analytic = |x: na::VectorView2<fx>| -> fx { (TAU * x[0]).cos() * (TAU * x[1]).cos() };
  let u_grad_analytic = |x: na::VectorView2<fx>| -> na::Vector2<fx> {
    na::Vector2::new(
      -TAU * (TAU * x[0]).sin() * (TAU * x[1]).cos(),
      -TAU * (TAU * x[0]).cos() * (TAU * x[1]).sin(),
    )
  };

  let galmat_laplacian = assemble_galmat(&mesh, elmat_laplacian);
  let galmat_mass = assemble_galmat(&mesh, elmat_mass);
  let galmat = galmat_laplacian + galmat_mass;
  let galvec = assemble_galvec(&mesh, |tria| elvec(tria, f));

  let mu = nas::factorization::CscCholesky::factor(&galmat)
    .unwrap()
    .solve(&galvec);

  let error_l2 = error_l2(&mesh, mu.as_view(), u_analytic);
  let error_h1s = error_h1s(&mesh, mu.as_view(), u_grad_analytic);

  dbg!(error_l2);
  dbg!(error_h1s);
}
