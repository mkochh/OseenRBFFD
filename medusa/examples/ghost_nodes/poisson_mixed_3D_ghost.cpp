#include <medusa/Medusa_fwd.hpp>
#include <medusa/bits/domains/GeneralFill.hpp>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>

using namespace mm;  // NOLINT
using namespace std;  // NOLINT
using namespace Eigen;  // NOLINT

/// Solving 3D Poisson's equation on unit ball
/// with mixed boundary conditions and ghost nodes for Neumann boundary.
/// http://e6.ijs.si/medusa/wiki/index.php/Ghost_nodes

using namespace mm;  // NOLINT

typedef Vec3d vec;  // Change to 2D or 1D

int main() {
    BallShape<vec> b1(0.0, 1.0);
    BallShape<vec> b2(1.0, 1.5);
    auto b = b1 - b2;
    double dx = 0.05;
    DomainDiscretization<vec> domain = b.discretizeBoundaryWithStep(dx);
    auto fn = [=](const vec&) { return dx; };
    GeneralFill<vec> fill; fill.seed(1337);
    fill(domain, fn);

    // Node indexes for interior, Dirichlet and Neumann BC.
    Range<int> interior = domain.interior();
    Range<int> neu, dir;
    for (int i : domain.boundary()) {
        double x = domain.pos(i, 0);
        if (x < 0) dir.push_back(i);
        else neu.push_back(i);
    }
    auto gh = domain.addGhostNodes(dx, 0, neu);

    // Find stencil.
    int N = domain.size();
    prn(N);
    int n = 25;
    domain.findSupport(FindClosest(n));

    // Declare approximation.
    std::cout << "Computing shapes..." << std::endl;
    Polyharmonic<double, 3> phs;
    RBFFD<decltype(phs), vec, ScaleToClosest> approx(phs, Monomials<vec>(2));
    auto storage = domain.computeShapes<sh::lap|sh::d1>(approx);

    Eigen::SparseMatrix<double, Eigen::RowMajor> M(N, N);
    M.reserve(Range<int>(N, n));
    Eigen::VectorXd rhs(N); rhs.setZero();

    auto op = storage.implicitOperators(M, rhs);

    std::cout << "Assembling..." << std::endl;
    for (int i : interior) {
        op.lap(i) = 1.0;  // Equation
    }
    for (int i : dir) {
        op.value(i) = 0.0;  // Dirichlet
    }
    for (int i : neu) {
        op.neumann(i, domain.normal(i)) = 0.0;  // Neumann
        // Equation holds also on the boundary, write it in the row corresponding to the ghost node.
        op.lap(i, gh[i]) = 1.0;
    }

    std::cout << "Solving..." << std::endl;
    // Solve the problem
    SparseLU<decltype(M)> solver;
    solver.compute(M);
    ScalarFieldd u = solver.solve(rhs);

    std::ofstream out_file("poisson_mixed_3D_ghost_data.m");
    out_file << "positions = " << domain.positions() << ";" << std::endl;
    out_file << "types = " << domain.types() << ";" << std::endl;
    out_file << "solution = " << u << ";" << std::endl;
    out_file.close();

    return 0;
}

