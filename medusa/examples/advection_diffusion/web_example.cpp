#include <medusa/Medusa_fwd.hpp>
#include <medusa/bits/domains/BasicRelax.hpp>
#include <medusa/bits/domains/GeneralFill.hpp>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

/// Example of advection-diffusion equation on http://e6.ijs.si/medusa/features

using namespace mm;  // NOLINT

template <int dim>
void solve() {
    typedef Vec<double, dim> vec;
    BallShape<vec> c(0.0, 1.0);
    BallShape<vec> hole(0.1, 0.3);
    BoxShape<vec> box({0.0, 0.0}, {1.0, 1.0});
    double dx = 0.001*dim*dim*dim;
    // DomainDiscretization<vec> domain = (c-hole).discretizeBoundaryWithStep(dx);
    DomainDiscretization<vec> domain = box.discretizeBoundaryWithStep(dx);
    GeneralFill<vec> fill;
    domain.fill(fill, dx);
    prn(domain.size());

    int N = domain.size();
    Monomials<vec> mon(2);
    domain.findSupport(FindClosest(12));

    std::cout << domain.support(domain.interior()[0]) << std::endl;

    RBFFD<Polyharmonic<double, 3>, vec, ScaleToClosest> approx({}, mon);
    auto storage = domain.template computeShapes<sh::lap|sh::d1>(approx);

    Eigen::SparseMatrix<double, Eigen::RowMajor> M(N, N);
    Eigen::VectorXd rhs(N); rhs.setZero();

    auto op = storage.implicitOperators(M, rhs);
    M.reserve(storage.supportSizes());

    for (int i : domain.interior()) {
        8.0*op.grad(i, -1) + 2.0*op.lap(i) = -1.0;
    }
    for (int i : domain.boundary()) {
        op.value(i) = 0.0;
    }

    Eigen::BiCGSTAB<decltype(M), Eigen::IncompleteLUT<double>> solver;
    solver.preconditioner().setDroptol(1e-4);
    solver.preconditioner().setFillfactor(20);
    solver.compute(M);
    ScalarFieldd u = solver.solve(rhs);

    std::ofstream out_file(format("web_example_%dd_data.m", dim));
    out_file << "positions = " << domain.positions() << ";" << std::endl;
    out_file << "solution = " << u << ";" << std::endl;
    out_file.close();
}


int main() {
    // solve<1>();
    solve<2>();
    // solve<3>();
    return 0;
}

