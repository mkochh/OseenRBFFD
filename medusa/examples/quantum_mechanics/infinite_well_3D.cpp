#include <medusa/Medusa_fwd.hpp>
#include <complex>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

/// We are solving Schrödinger equation for a 3D infinite potential well in a unit cube.
/// http://e6.ijs.si/medusa/wiki/index.php/Schr%C3%B6dinger_equation#3D_infinite_cubic_well

using namespace mm;  // NOLINT
using namespace std;  // NOLINT
using namespace Eigen;  // NOLINT

constexpr std::complex<double> operator""_i(long double d) {
    return std::complex<double>{0.0, static_cast<double>(d)};
}

int main() {
    // Set parameters.
    double step = 0.05;
    int n = 9;  // support size
    double T = 0.1;
    double dt = 1e-5;

    // Create the domain and discretize it.
    BoxShape<Vec3d> box(0.0, 1.0);
    DomainDiscretization<Vec3d> domain = box.discretizeWithStep(step);

    // Find support for the nodes.
    domain.findSupport(FindClosest(n));
    int N = domain.size();

    // Construct the approximation engine.
    WLS<Gaussians<Vec3d>, NoWeight<Vec3d>, ScaleToFarthest> wls({9, 30});

    // Prepare operators and matrix.
    auto storage = domain.computeShapes<sh::lap>(wls);  // shapes are computed using our wls

    SparseMatrix<complex<double>, RowMajor> M(N, N);
    M.reserve(storage.supportSizes());
    VectorXcd rhs = VectorXd::Zero(N);
    VectorXcd psi(N);
    auto op = storage.implicitOperators(M, rhs);


    // Set initial state.
    for (int i : domain.interior()) {
        double x = domain.pos(i, 0);
        double y = domain.pos(i, 1);
        double z = domain.pos(i, 2);
        rhs(i) = sin(PI * x) * sin(PI * y) * sin(PI * z);
    }

    // Set equation on interior.
    for (int i : domain.interior()) {
        op.value(i) + (-1.0_i) * dt * op.lap(i) = rhs(i);
    }

    // Set equation on boundary.
    for (int i : domain.boundary()) {
        op.value(i) = 0;
    }

    BiCGSTAB<decltype(M), IncompleteLUT<complex<double>>> solver;
    solver.compute(M);

    // Time stepping.
    int steps = iceil(T/dt);
    for (int t = 1; t < steps; ++t) {
        // solve matrix system
        psi = solver.solve(rhs);
        // update rhs
        rhs = psi;
    }

    // Write the solution into file.
    HDF hdf("infinite_well_3D.h5", HDF::DESTROY);
    hdf.writeDouble2DArray("pos", domain.positions());
    hdf.writeDoubleArray("rsol", psi.real());
    hdf.writeDoubleArray("csol", psi.imag());
    hdf.close();

    return 0;
}
