#include <medusa/Medusa_fwd.hpp>
#include <medusa/bits/domains/BasicRelax.hpp>
#include <medusa/bits/domains/GeneralFill.hpp>
#include "Eigen/SparseCore"
#include "Eigen/IterativeLinearSolvers"

using namespace mm;  // NOLINT

int main() {
    // Construct a box [0,1]*[0,1]
    BoxShape<Vec2d> box(0.0, 1.0);
    BallShape<Vec2d> c0({0.25, 0.25}, 0.05);
    BallShape<Vec2d> c1({0.25, 0.75}, 0.05);
    BallShape<Vec2d> c2({0.75, 0.75}, 0.05);
    BallShape<Vec2d> c3({0.75, 0.25}, 0.05);

    // Discretize the domain
    double step = 0.01;
    DomainDiscretization<Vec2d> domain = (box-c0 - c1 - c2 - c3).discretizeBoundaryWithStep(step);

    // Fill the domain
    GeneralFill<Vec2d> fill;
    domain.fill(fill, step);

    // Relax the domain
    BasicRelax relax;
    relax.iterations(10).initialHeat(0.8).numNeighbours(5)
            .projectionType(BasicRelax::DO_NOT_PROJECT);
    relax(domain, step);

    // Find support nodes
    FindClosest find_support(9);
    domain.findSupport(find_support);

    // Define approximation engine
    WLS<Gaussians<Vec2d>, NoWeight<Vec2d>, ScaleToFarthest,
            Eigen::LLT<Eigen::MatrixXd>> wls({9, 30.0});

    // Compute shapes
    auto storage = domain.computeShapes<sh::lap>(wls);

    // Construct Sparse Matrix to store the shapes
    int N = domain.size();
    Eigen::SparseMatrix<double, Eigen::RowMajor> M(N, N);
    M.reserve(storage.supportSizes());

    // Construct vector to store the right hand side of the equation
    Eigen::VectorXd rhs(N); rhs.setZero();
    auto op = storage.implicitOperators(M, rhs);

    for (int i : domain.interior()) {
        op.lap(i) = 0;
    }
    for (int i : (domain.types() == -1)) {
        double y = domain.pos(i, 1);
        double x = domain.pos(i, 0);

        if (x > 0.5) {
            if (y > 0.5) {
                op.value(i) = 1;
            } else {
                op.value(i) = -1;
            }
        } else {
            if (y > 0.5) {
                op.value(i) = -1;
            } else {
                op.value(i) = 1;
            }
        }

        if (x == 0.0) {
            op.value(i) = 0;
        }
    }
    for (int i : (domain.types() == -2)) {
        op.value(i) = 0;
    }
    for (int i : (domain.types() == -3)) {
        op.value(i) = 0;
    }
    for (int i : (domain.types() == -4)) {
        op.value(i) = 0;
    }

    Eigen::BiCGSTAB<decltype(M), Eigen::IncompleteLUT<double>> solver;
    solver.compute(M);
    ScalarFieldd u = solver.solve(rhs);

    HDF hdf("quadrupole.h5", HDF::DESTROY);
    hdf.writeDomain("domain", domain);
    hdf.writeDoubleArray("sol", u);
    hdf.close();

    return 0;
}
