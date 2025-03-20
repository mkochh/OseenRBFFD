#ifndef RBFFD_DISCRETIZATION_HEADER
#define RBFFD_DISCRETIZATION_HEADER

// #include "../../medusa/include/medusa/Medusa.hpp"
#include <medusa/Medusa_fwd.hpp>
#include "../../medusa/include/medusa/bits/domains/PolyhedronShape.hpp"
#include "../auxiliaries/polyhedron_integration.h"

namespace mm {

class OseenDiscretizationBetter {
  
  public:
    DomainDiscretization<Vec3d> d_u;
    DomainDiscretization<Vec3d> d_p;
    DomainDiscretization<Vec3d> d_grad;
    DomainDiscretization<Vec3d> d_div;
    DomainDiscretization<Vec3d> d_conv;
    DomainDiscretization<Vec3d> d_u_int;
    Range<int> idxs_p;
    Range<int> idxs_u;
    Range<int> idxs_ui;
    int N_u;
    int N_ui;
    int N_ub;
    int N_p;
    int N_dofs;
    int poly_lap;
    int poly_conv;
    int poly_grad;
    int poly_div;
    double nu;
    Eigen::VectorXd constraint;
    Eigen::VectorXd exact_solution;
    mm::Range<mm::Vec3d> conv_vec3d;
    mm::Range<mm::Vec3d> dir_bnd_vec3d;
    mm::Range<mm::Vec3d> rhs_vec3d;
    int seed;
    int n[6] = {0, 10, 22, 43, 74, 116}; // stencil sizes for different deg of polynomial augmentation
    int l[6] = {0, 1, 2, 3, 4, 5}; // degree of polynomial augmentation
    int k[6] = {0, 1, 3, 3, 3, 5}; // degree of polyharmonic spline for different deg of polynomial augmentation

    OseenDiscretizationBetter(PolyhedronShape<Vec3d> shape = PolyhedronShape<Vec3d>::fromOFF("OFF_Files/cube.off"),
                       double dx_u = 0.1, double dx_p = 0.2, 
                       int poly_grad = 3, int poly_conv = 3, int poly_lap = 4, int poly_div = 3, double nu = 1e-2, int seed = 0)
                       : d_u(shape)
                       , d_p(shape)
                       , d_grad(shape)
                       , d_div(shape)
                       , d_conv(shape)
                       , d_u_int(shape)
                       , poly_lap(poly_lap)
                       , poly_conv(poly_conv)
                       , poly_grad(poly_grad)
                       , poly_div(poly_div)
                       , nu(nu)
                       , seed(seed)
                      { initialize(shape, dx_u, dx_p); }

    void determineSupports() {
      FindClosest supp_lap(n[poly_lap]);
      supp_lap.forNodes(idxs_ui); // don't need stencils at boundary
      d_u.findSupport(supp_lap); // find support for lap within d_u
      FindClosest supp_conv(n[poly_conv]);
      supp_conv.forNodes(idxs_ui); // don't need stencils at boundary
      d_conv.findSupport(supp_conv); // find support for conv within d_conv(=d_u)
      FindClosest supp_grad(n[poly_grad]); // find support for grad over nodes in d_p at points in d_u
      supp_grad.forNodes(idxs_ui).searchAmong(idxs_p).forceSelf(false);
      d_grad.findSupport(supp_grad);
      FindClosest supp_divu(n[poly_div]); // find support for div u over nodes in d_u at points in d_p
      supp_divu.forNodes(idxs_p).searchAmong(idxs_u).forceSelf(false);
      d_div.findSupport(supp_divu);
    }

    void setConstraint(std::string OFF_file, int max_order = 12) {
      constraint = compute_pressure_constraint(OFF_file, d_p, max_order);
    }

    void setConvection(int which_convection = 0);

    void setDirBoundary(int which_dir_bnd = 1);

    void setRHS(int which_rhs = 0);

    void computeExactSolution(int which_solution, std::string domain_name);

    void setSolution(int which_solution, std::string domain_name);

  protected:
    void initialize(PolyhedronShape<Vec3d> shape, double dx_u, double dx_p);
};

void OseenDiscretizationBetter::initialize(PolyhedronShape<Vec3d> shape, double dx_u, double dx_p) {
  // create point set for Laplacian of velocity
  d_u = shape.discretizeBoundaryWithStep(dx_u);
  GeneralFill<Vec3d> fill;
  fill.seed(seed);
  d_u.fill(fill,dx_u); // generate nodes

  // laplacian and velocity use same point set, but separate them to use different stencils
  d_conv = d_u;

  // create point set for pressure
  d_p = shape.discretizeBoundaryWithStep(dx_p*0.5); // discretize boundary, this is removed later and only needed to place interior and ghost points
  Range<int> bd = d_p.boundary(); // save initial boundary node indices
  // dx_p/2.5 can be changed to place ghost nodes closer or further away from boundary
  // respectively pressure nodes will be further or closer to the boundary
  Range<int> gh = d_p.addGhostNodes(dx_p/2.5,-2); // add ghost nodes as boundary nodes
  d_p.removeNodes(bd); // remove initial boundary nodes
  d_p.fill(fill,dx_p); // fill domain defined by new boundary nodes
  d_p.removeBoundaryNodes(); // remove new boundary nodes
  // comment out 349-356 and uncomment 358-359 to have pressure nodes equal interior velocity nodes
  // d_p = d_u;
  // d_p.removeBoundaryNodes();

  // create combined point sets for gradient and divergence and a set of interior velocity nodes
  d_grad = d_u;
  idxs_p = d_grad.addNodes(d_p);
  idxs_u = d_u.all();
  idxs_ui = d_u.interior();
  d_div = d_grad;
  d_u_int = d_u;
  d_u_int.removeBoundaryNodes();

  N_u = d_u.size();
  N_ui = d_u.interior().size();
  N_ub = d_u.boundary().size();
  N_p = d_p.size();
  N_dofs = 3*N_ui+N_p+1;

  mm::Range<mm::Vec3d> temp(N_u,0);
  conv_vec3d = temp;
  dir_bnd_vec3d = temp;
  rhs_vec3d = temp;

  constraint = Eigen::VectorXd::Ones(N_p);
  exact_solution = Eigen::VectorXd::Zero(N_dofs);
}

void OseenDiscretizationBetter::setConvection(int which_convection) {
  switch (which_convection) {
    // simple convection in x-direction
    case 0:
      for (int i : idxs_ui)
        conv_vec3d[i] = {1.0, 0.0, 0.0};
      break;
    // more complicated convection
    case 1: {
      double scaling = 6.40936351829;
      for (int i : idxs_ui) {
        double x = d_u.pos(i,0);
        double y = d_u.pos(i,1);
        double z = d_u.pos(i,2);
        conv_vec3d[i] = {scaling*2.0*x*(1.0-x)*(2.0*y-1)*z, -scaling*(2.0*x-1)*y*(1.0-y), -scaling*(2.0*x-1)*(2.0*y-1.0)*z*(1.0-z)};
      }
      }
      break;
    // no convection (Stokes equations)
    case 2:
      for (int i : idxs_ui)
        conv_vec3d[i] = {0.0, 0.0, 0.0};
      break;
      
    default:
      break;
  }
}

void OseenDiscretizationBetter::setDirBoundary(int which_dir_bnd) {
  switch (which_dir_bnd) {
    case 0:
      for (int i : d_u.boundary())
        dir_bnd_vec3d[i] = {0.0, 0.0, 0.0};
      break;
    case 1:
      for (int i : d_u.boundary())
        dir_bnd_vec3d[i] = {1.0, 0.0, 0.0};
      break;
    case 2:
      for (int i : d_u.boundary())
        dir_bnd_vec3d[i] = {1.0, 1.0, 1.0};
      break;
    case 3:
      for (int i : d_u.boundary()) {
        double x = d_u.pos(i,0);
        double y = d_u.pos(i,1);
        double z = d_u.pos(i,2);
        dir_bnd_vec3d[i] = {cos(x)*sin(y), x*cos(y), z*(x+sin(x))*sin(y)};
      }
      break;
      
    default:
      break;
  }
}

// set the right hand side according to the chosen exact solution
void OseenDiscretizationBetter::setRHS(int which_rhs) {
  switch (which_rhs) {
  case 0:
      for (int i : idxs_ui) {
          double x = d_u.pos(i,0);
          double y = d_u.pos(i,1);
          rhs_vec3d[i] = {PI*cos(PI*x)*cos(PI*y),
                            -PI*sin(PI*x)*sin(PI*y),
                            0.0};
      }
      break;
  case 1:
      for (int i : idxs_ui) {
          double x = d_u.pos(i,0);
          double y = d_u.pos(i,1);
          double z = d_u.pos(i,2);
          rhs_vec3d[i] = {2.0*nu*cos(x)*sin(y)-conv_vec3d[i](0)*sin(x)*sin(y)+conv_vec3d[i](1)*cos(x)*cos(y) + PI*cos(PI*x)*cos(PI*y),
                            nu*x*cos(y)+conv_vec3d[i](0)*cos(y)-conv_vec3d[i](1)*x*sin(y) - PI*sin(PI*x)*sin(PI*y),
                            nu*z*(x+2.0*sin(x))*sin(y)+conv_vec3d[i](0)*z*(cos(x)+1)*sin(y)+conv_vec3d[i](1)*z*(x+sin(x))*cos(y)+conv_vec3d[i](2)*(x+sin(x))*sin(y)};
      }
      break;

  case 2:
      for (int i : idxs_ui) {
          double x = d_u.pos(i,0);
          double y = d_u.pos(i,1);
          double z = d_u.pos(i,2);
          rhs_vec3d[i] = {(y-0.5)*(z-0.5), (x-0.5)*(z-0.5), (x-0.5)*(y-0.5)};
      }
      break;

  case 3:
      for (int i : idxs_ui) {
          double x = d_u.pos(i,0);
          double y = d_u.pos(i,1);
          double z = d_u.pos(i,2);
          rhs_vec3d[i] = {2.0*nu*cos(x)*sin(y)-conv_vec3d[i](0)*sin(x)*sin(y)+conv_vec3d[i](1)*cos(x)*cos(y) + (y-0.5)*(z-0.5),
                            nu*x*cos(y)+conv_vec3d[i](0)*cos(y)-conv_vec3d[i](1)*x*sin(y) + (x-0.5)*(z-0.5),
                            nu*z*(x+2.0*sin(x))*sin(y)+conv_vec3d[i](0)*z*(cos(x)+1)*sin(y)+conv_vec3d[i](1)*z*(x+sin(x))*cos(y)+conv_vec3d[i](2)*(x+sin(x))*sin(y) + (x-0.5)*(y-0.5)};
      }
      break;
  
  default:
      break;
  }
}

// set the exact solution
void OseenDiscretizationBetter::computeExactSolution(int which_solution, std::string domain_name) {
  switch (which_solution)
  {
  case 0:
    for(int i : idxs_ui) {
        exact_solution(i-N_ub) = 0.0;
        exact_solution(i+N_ui-N_ub) = 0.0;
        exact_solution(i+2*N_ui-N_ub) = 0.0;
    }
    for(int i : d_p.all()) {
        double x = d_p.pos(i,0);
        double y = d_p.pos(i,1);
        if (domain_name == "bunny") {
            exact_solution(i+3*N_ui) = std::sin(PI*x)*std::cos(PI*y) - 0.260568484878944; // integral of this over bunny domain is zero
        } else {
            exact_solution(i+3*N_ui) = std::sin(PI*x)*std::cos(PI*y); // integral of this over cube and narrowing domain is zero
        }
    }
    break;
  case 1:
    for(int i : idxs_ui) {
        exact_solution(i-N_ub) = 1.0;
        exact_solution(i+N_ui-N_ub) = 0.0;
        exact_solution(i+2*N_ui-N_ub) = 0.0;
    }
    for(int i : d_p.all()) {
        double x = d_p.pos(i,0);
        double y = d_p.pos(i,1);
        if (domain_name == "bunny") {
            exact_solution(i+3*N_ui) = std::sin(PI*x)*std::cos(PI*y) - 0.260568484878944; // integral of this over bunny domain is zero
        } else {
            exact_solution(i+3*N_ui) = std::sin(PI*x)*std::cos(PI*y); // integral of this over cube and narrowing domain is zero
        }
    }
    break;
  case 2:
    for(int i : idxs_ui) {
        exact_solution(i-N_ub) = 1.0;
        exact_solution(i+N_ui-N_ub) = 1.0;
        exact_solution(i+2*N_ui-N_ub) = 1.0;
    }
    for(int i : d_p.all()) {
        double x = d_p.pos(i,0);
        double y = d_p.pos(i,1);
        if (domain_name == "bunny") {
            exact_solution(i+3*N_ui) = std::sin(PI*x)*std::cos(PI*y) - 0.260568484878944; // integral of this over bunny domain is zero
        } else {
            exact_solution(i+3*N_ui) = std::sin(PI*x)*std::cos(PI*y); // integral of this over cube and narrowing domain is zero
        }
    }
    break;
  case 3:
    for(int i : idxs_ui) {
        double x = d_u.pos(i,0);
        double y = d_u.pos(i,1);
        double z = d_u.pos(i,2);
        exact_solution(i-N_ub) = cos(x)*sin(y);
        exact_solution(i+N_ui-N_ub) = x*cos(y);
        exact_solution(i+2*N_ui-N_ub) = z*(x+sin(x))*sin(y);
    }
    for(int i : d_p.all()) {
        double x = d_p.pos(i,0);
        double y = d_p.pos(i,1);
        if (domain_name == "bunny") {
            exact_solution(i+3*N_ui) = std::sin(PI*x)*std::cos(PI*y) - 0.260568484878944; // integral of this over bunny domain is zero
        } else {
            exact_solution(i+3*N_ui) = std::sin(PI*x)*std::cos(PI*y); // integral of this over cube and narrowing domain is zero
        }
    }
    break;

  case 4:
    for(int i : idxs_ui) {
        exact_solution(i-N_ub) = 1.0;
        exact_solution(i+N_ui-N_ub) = 1.0;
        exact_solution(i+2*N_ui-N_ub) = 1.0;
    }
    for(int i : d_p.all()) {
        double x = d_p.pos(i,0);
        double y = d_p.pos(i,1);
        double z = d_p.pos(i,2);
        if (domain_name == "bunny") {
            exact_solution(i+3*N_ui) = (x-0.5)*(y-0.5)*(z-0.5) - 0.001363726470335; // integral of this over cube and narrowing domain is zero
        } else {
            exact_solution(i+3*N_ui) = (x-0.5)*(y-0.5)*(z-0.5); // integral of this over cube and narrowing domain is zero
        }
    }
    break;

  case 5:
    for(int i : idxs_ui) {
        double x = d_u.pos(i,0);
        double y = d_u.pos(i,1);
        double z = d_u.pos(i,2);
        exact_solution(i-N_ub) = cos(x)*sin(y);
        exact_solution(i+N_ui-N_ub) = x*cos(y);
        exact_solution(i+2*N_ui-N_ub) = z*(x+sin(x))*sin(y);
    }
    for(int i : d_p.all()) {
        double x = d_p.pos(i,0);
        double y = d_p.pos(i,1);
        double z = d_p.pos(i,2);
        if (domain_name == "bunny") {
            exact_solution(i+3*N_ui) = (x-0.5)*(y-0.5)*(z-0.5) - 0.001363726470335; // integral of this over cube and narrowing domain is zero
        } else {
            exact_solution(i+3*N_ui) = (x-0.5)*(y-0.5)*(z-0.5); // integral of this over cube and narrowing domain is zero
        }
    }
    break;
  
  default:
    break;
  }
}

// sets Dirichlet boundary condition, RHS and exact solution such that they fit together
void OseenDiscretizationBetter::setSolution(int which_solution, std::string domain_name) {
  switch (which_solution)
  {
  case 0:
    setDirBoundary(0);
    setRHS(0);
    computeExactSolution(0, domain_name);
    break;
  case 1:
    setDirBoundary(1);
    setRHS(0);
    computeExactSolution(1, domain_name);
    break;
  case 2:
    setDirBoundary(2);
    setRHS(0);
    computeExactSolution(2, domain_name);
    break;
  case 3:
    setDirBoundary(3);
    setRHS(1);
    computeExactSolution(3, domain_name);
    break;
  case 4:
    setDirBoundary(2);
    setRHS(2);
    computeExactSolution(4, domain_name);
    break;
  case 5:
    setDirBoundary(3);
    setRHS(3);
    computeExactSolution(5, domain_name);
    break;
  
  default:
    break;
  }
}

} // namespace mm

#endif