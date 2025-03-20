#include <medusa/Medusa.hpp>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

/// Basic medusa example, we are solving 2D Poisson's equation on a NURBS domain
/// with Dirichlet boundary conditions.
/// http://e6.ijs.si/medusa/wiki/index.php/NURBS_domains

using namespace mm;  // NOLINT

int main() {
    // Points to be used in NURBS creation.
    Range<Vec2d> pts{{210, -132.5}, {205, -115}, {125, -35}, {85, -61}, {85, -61}, {85, -61},
                     {80, -65}, {75, -58}, {75, -58}, {65, 45}, {25, 25}, {-15, -16}, {-15, -16},
                     {-15, -16}, {-35, -28}, {-40, -32.375}, {-40, -32.375}, {-43, -35}, {-27, -40},
                     {-5, -40}, {-5, -40}, {50, -38}, {35, -65}, {15, -75}, {15, -75}, {-15, -95},
                     {-20, -140}, {45, -146}, {45, -146}, {95, -150}, {215, -150}, {210, -132.5}};

    // Calculate NURBS patches' control points, weights, knot vector and order.
    int p = 3;
    Range<double> weights{1, 1, 1, 1};
    Range<double> knot_vector{0, 0, 0, 0, 1, 1, 1, 1};
    Range<NURBSPatch<Vec2d, Vec1d>> patches;

    for (int i = 0; i < 8; i++) {
        Range<Vec2d> control_points;
        for (int j = 0; j < 4; j++) {
            control_points.push_back(pts[4 * i + j]);
        }

        NURBSPatch<Vec2d, Vec1d> patch(control_points, weights, {knot_vector}, {p});
        patches.push_back(patch);
    }

    // Create NURBS shape from a range of patches.
    NURBSShape<Vec2d, Vec1d> shape(patches);

    // Fill the domain.
    double h = 0.5;
    DomainDiscretization<Vec2d> domain = shape.discretizeBoundaryWithStep(h);

    GeneralFill<Vec2d> gf;
    gf(domain, h);

    // Construct the approximation engine.
    int m = 2;  // basis order
    Monomials<Vec2d> mon(m);
    RBFFD<Polyharmonic<double, 3>, Vec2d> approx({}, mon);

    // Find support for the nodes.
    int N = domain.size();
    domain.findSupport(FindClosest(mon.size()));

    // Compute the shapes (we only need the Laplacian).
    auto storage = domain.computeShapes<sh::lap>(approx);

    Eigen::SparseMatrix<double, Eigen::RowMajor> M(N, N);
    Eigen::VectorXd rhs(N);
    rhs.setZero();
    M.reserve(storage.supportSizes());

    // Construct implicit operators over our storage.
    auto op = storage.implicitOperators(M, rhs);

    for (int i : domain.interior()) {
        // set the case for nodes in the domain
        op.lap(i) = 20.0;
    }
    for (int i : domain.boundary()) {
        // enforce the boundary conditions
        op.value(i) = 0.0;
    }

    Eigen::BiCGSTAB<decltype(M), Eigen::IncompleteLUT<double>> solver;
    solver.compute(M);
    ScalarFieldd u = solver.solve(rhs);

    // Write the solution into file.
    std::ofstream out_file("poisson_NURBS_2D_data.m");
    out_file << "positions = " << domain.positions() << ";" << std::endl;
    out_file << "solution = " << u << ";" << std::endl;
    out_file.close();

    return 0;
}
