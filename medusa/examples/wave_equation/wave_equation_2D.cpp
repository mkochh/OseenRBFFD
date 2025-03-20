#include <medusa/Medusa_fwd.hpp>
#include <medusa/bits/domains/GeneralFill.hpp>
#include <Eigen/Sparse>
#include <cmath>

/// 2D Wave equation problem on bounded domain.
/// http://e6.ijs.si/medusa/wiki/index.php/Wave_equation

using namespace mm;  // NOLINT
using namespace Eigen;  // NOLINT

// Helper function calculating linear interpolation/extrapolation y(x) based on
// two data points (x1, y1), (x2, y2) at point x.
double linear(double x1, double y1, double x2, double y2, double x) {
    double k = (y2 - y1) / (x2 - x1);
    double n = y1 -k*x1;
    return k*x + n;
}

int main() {
    // parameters
    double inner_radius = 0.05;
    double outer_radius = 1.0;
    int density = 40;
    double A = 0.2;
    int n = 12;  // support size
    int m = 3;  // monomial basis of second order, i.e. 3 monomials
    int fill_seed = 123;
    double sigma = 1;

    double omega = 128;  // frequency of the source
    double v = 8;  // wave velocity
    double time = 0.5;  // length of calculation
    int t_factor = 500;

    double dt = 1e-6;  // time step
    int t_steps = static_cast<int> (time / dt);
    double dx = outer_radius / density;

    // Lambda function for setting the density of nodes
    auto fill_density = [=](const Vec2d& p) {
        double r = p.norm();
        double default_value = dx;
        double dens = default_value;
        double r1 = 15*inner_radius;
        double r2 = 0.8*outer_radius;
        if (r < r1) dens = linear(inner_radius, 0.8*default_value, r1, default_value, r );
        if (r > r2) dens = linear(r2, default_value, outer_radius, 0.8* default_value, r);
        return dens;
    };

    // Create output file
    HDF hdf_out("wave_equation_2D.h5", HDF::DESTROY);

    // extra boundary labels
    int CENTRE = -10;

    // Prepare domain

    // build circle domain
    BallShape<Vec2d> domain({0, 0}, outer_radius);
    auto discretization = domain.discretizeBoundaryWithStep(dx);

    // build source
    BallShape<Vec2d> empty({0, 0}, inner_radius);
    auto discretization_empty = empty.discretizeBoundaryWithStep(dx, CENTRE);

    // substract the source domain
    discretization -= discretization_empty;


    GeneralFill<Vec2d> fill;
    fill.seed(fill_seed);
    discretization.fill(fill, fill_density);

    // find support
    FindClosest find_support(n);
    discretization.findSupport(find_support);

    int domain_size = discretization.size();

    Range<int> interior = discretization.types() > 0;
    Range<int> boundary = discretization.types() < 0;

    // Prepare operators and matrix
    WLS<Monomials<Vec2d>, GaussianWeight<Vec2d>, ScaleToClosest>
            approx(m - 1 , sigma);
    auto storage = discretization.computeShapes(approx);  // Shape functions are computed.

    SparseMatrix<double> M(domain_size, domain_size);
    BiCGSTAB<SparseMatrix<double, RowMajor>> solver;
    M.reserve(Range<int>(domain_size, n));
    VectorXd rhs = VectorXd::Zero(domain_size);  // set empty vector for rhs

    auto op = storage.implicitOperators(M, rhs);  // All nodes, including boundary

    VectorXd E0 = VectorXd::Zero(domain_size);  // 0-th step
    VectorXd E1 = VectorXd::Zero(domain_size);  // 1-st step

    // Set equation on interior
    for (int i : interior) {
        -(v*v*dt*dt) * op.lap(i) + op.value(i)= 2*E1(i)-E0(i);  // laplace in interior
    }
    // Set boundary conditions
    for (int i : boundary) {
        if (discretization.types()[i] == CENTRE) {
            op.value(i) = A * std::sin(omega*dt*2);
        } else {
            op.value(i) = 0.0;
        }
    }
    hdf_out.writeSparseMatrix("M", M);
    M.makeCompressed();
    solver.compute(M);


    // Time stepping
    int tt;
    int t_save = 0;
    hdf_out.writeDouble2DArray("pos", discretization.positions());
    hdf_out.openGroup("/step" + std::to_string(t_save));
    hdf_out.writeDoubleAttribute("TimeStep", t_save);
    hdf_out.writeDoubleAttribute("time", 0.0);
    hdf_out.writeDoubleArray("E", E1);
    ++t_save;

    for (tt = 2; tt < t_steps; ++tt) {
        // Solve matrix system
        VectorXd E2 = solver.solve(rhs);

        // Update previous states
        E0 = E1;
        E1 = E2;

        if (tt%(t_factor) == 0) {
            // Saving state
            hdf_out.openGroup("/step" + std::to_string(t_save));
            hdf_out.writeDoubleAttribute("TimeStep", t_save);
            std::cout<< "Time step " << t_save << " of " << t_steps/t_factor << "." << std::endl;
            ++t_save;
            hdf_out.writeDoubleArray("E", E2);
            hdf_out.writeDoubleAttribute("time", dt * tt);
        }

        // Update rhs on interior
        for (int i : interior) {
            rhs(i) = 2*E1(i)-E0(i);
        }

        // Update rhs on boundary
        for (int i : (discretization.types() == CENTRE)) {
            rhs(i) = A * std::sin(omega*dt*tt);
        }
    }
    hdf_out.closeFile();
}
