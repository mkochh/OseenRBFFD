#include <medusa/Medusa_fwd.hpp>
#include <medusa/bits/domains/GeneralFill.hpp>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>
#include <cmath>

// Basic coupled domains example 2D.
// http://e6.ijs.si/medusa/wiki/index.php/Coupled_domains

using namespace mm;  // NOLINT

int main() {
    double r1 = 0.5;  // inner radius
    double r2 = 1;  // outer radius
    double lam1 = 5;  // outer lambda
    double lam2 = 1;  // inner lambda
    double h =  0.01;  // nodal spacing
    double q0 = 1;  // volumetric heat source in interior
    int n = 9;  // support size

    BallShape<Vec2d> inner_shape(0, r1);
    BallShape<Vec2d> outer_circle(0, r2);
    auto outer_shape = outer_circle - inner_shape;

    DomainDiscretization<Vec2d> outer = outer_shape.discretizeBoundaryWithStep(h);
    DomainDiscretization<Vec2d> inner(inner_shape);

    // Indices of nodes in outer domain that constitute outer and common boundary.
    auto common_bnd = outer.positions().filter([=](const Vec2d& p) {
                            return std::abs(p.norm() - r1) < 1e-5; });

    auto outer_bnd = outer.positions().filter([=](const Vec2d& p) {
                            return std::abs(p.norm() - r2) < 1e-5; });

    // Maps node indices of outer domain to their corresponding indices in the inner domain.
    Range<int> outer_to_inner(outer.size(), -1);
    for (int i : common_bnd) {
        int idx = inner.addBoundaryNode(outer.pos(i), -2, -outer.normal(i));
        outer_to_inner[i] = idx;
    }

    GeneralFill<Vec2d> fill_engine; fill_engine.seed(1);
    outer.fill(fill_engine, h);  // this is the annulus
    inner.fill(fill_engine, h);  // this is the inner circle

    HDF out("poisson_coupled_domains.h5", HDF::DESTROY);
    out.writeIntArray("map", outer_to_inner);
    out.writeDomain("inner", inner);
    out.writeDomain("outer", outer);
    out.writeDoubleAttribute("lam1", lam1);
    out.writeDoubleAttribute("lam2", lam2);
    out.close();


    // Find support
    inner.findSupport(FindClosest(n));
    outer.findSupport(FindClosest(n));

    ////////////////////////////// SOLVING ///////////////////////////////////

    WLS<Gaussians<Vec2d>, GaussianWeight<Vec2d>, ScaleToFarthest> wls({9, 30}, 1.0);

    int N_inner = inner.size();
    int N_outer = outer.size();

    Eigen::SparseMatrix<double, Eigen::RowMajor> M(N_inner+N_outer, N_inner+N_outer);
    Eigen::VectorXd rhs(N_inner+N_outer);
    rhs.setZero();

    auto storage_inner = inner.computeShapes<sh::lap|sh::d1>(wls);
    auto storage_outer = outer.computeShapes<sh::lap|sh::d1>(wls);

    auto op_inner = storage_inner.implicitOperators(M, rhs);
    auto op_outer = storage_outer.implicitOperators(M, rhs);
    op_outer.setRowOffset(N_inner);
    op_outer.setColOffset(N_inner);

    Range<int> reserve = storage_inner.supportSizes() + storage_outer.supportSizes();
    // Compute row sizes on the common boundary
    for (int i : common_bnd) {
        reserve[outer_to_inner[i]] = 2;  // equality of values needs 2 slots
        reserve[N_inner+i] = storage_inner.supportSize(outer_to_inner[i])
                + storage_outer.supportSize(i);  // equality of fluxes needs two shapes
    }
    M.reserve(reserve);

    ///  CASE
    for (int i : inner.interior()) {
        -lam2*op_inner.lap(i) = q0;  // laplace in interior of inner domain
    }
    for (int i : outer.interior()) {
        lam1*op_outer.lap(i) = 0;  // laplace in interior of outer domain
    }
    for (int i : outer_bnd) {
        op_outer.value(i) = 0.0;  // outer Dirichlet BC
    }
    for (int i : common_bnd) {
        int j = outer_to_inner[i];
        // continuity over shared boundary
        op_inner.value(j) + -1.0*op_outer.value(i, -N_inner+j) = 0.0;
        lam2*op_inner.neumann(j, inner.normal(j), N_inner+i) +
        lam1*op_outer.neumann(i, outer.normal(i)) = 0.0;  // flux
    }

    // reserve and construct RHS
    out.atomic().writeSparseMatrix("M", M);
    out.atomic().writeDoubleArray("rhs", rhs);

    Eigen::BiCGSTAB<decltype(M), Eigen::IncompleteLUT<double>> solver;
    solver.preconditioner().setFillfactor(50);
    solver.preconditioner().setDroptol(1e-5);
    solver.setMaxIterations(100);
    solver.setTolerance(1e-15);
    solver.compute(M);
    Eigen::VectorXd sol = solver.solve(rhs);
    assert_msg(solver.info() == Eigen::Success, "Coupled domain example did not converge.");
    prn(solver.iterations())
    prn(solver.error())

    out.atomic().writeDoubleArray("sol", sol);

    return 0;
}
