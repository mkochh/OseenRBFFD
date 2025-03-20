#include <medusa/Medusa_fwd.hpp>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>


/// Basic medusa example, we are solving 2D lid driven case on unit square
/// with Dirichlet boundary conditions.
/// http://www-e6.ijs.si/medusa/wiki/index.php/Lid_driven_cavity

using namespace mm;  // NOLINT

int main() {
    // case setup
    double dt = 0.001;
    double time = 0.5;
    int n = 9;
    double sigma_b = 20;
    double dx = 0.02;
    double mu = 0.01;
    double rho = 1;
    double dl = 2;
    // Create the domain and discretize it
    BoxShape<Vec2d> box(0.0, 1.0);
    DomainDiscretization<Vec2d> domain = box.discretizeWithStep(dx);
    // boundary maps
    Range<int> top = domain.types() == -4;
    // fields declaration
    VectorField2d u_1(domain.size());
    ScalarFieldd p = ScalarFieldd::Zero(domain.size());
    u_1.setZero();
    // dirichlet boundary conditions - lid is set to {1,0}, other remain zero
    u_1[top] = {1, 0};
    VectorField2d u_2 = u_1;
    // Find support for the nodes
    domain.findSupport(FindClosest(n));  // the support for each node is the closest 9 nodes

    // Construct the approximation engine, in this case a weighted least squares using monomials as
    // basis functions, no weight, and scale to farthest
    WLS<Gaussians<Vec2d>, NoWeight<Vec2d>,
            ScaleToFarthest> wls({n, sigma_b});

    // compute the shapes (we only need the Laplacian) using our WLS
    auto storage = domain.computeShapes<sh::lap | sh::grad>(wls);
    // construct explicit operators
    auto op_v = storage.explicitVectorOperators();
    auto op_s = storage.explicitOperators();

    // time loop
    for (int time_step = 0; time_step < std::floor(time / dt); ++time_step) {
        // Navier-Stokes
        for (int i : domain.interior()) {
            u_2[i] = u_1[i] + dt * (-1 / rho * op_s.grad(p, i) + mu / rho * op_v.lap(u_1, i)
                                    - op_v.grad(u_1, i) * u_1[i]);
        }
        for (int i : domain.all()) {
            p[i] = p[i] - dl * dl * dt * rho * op_v.div(u_2, i);
        }
        u_1.swap(u_2);
    }

    // Write the solution to file
    std::ofstream out_file("lid_driven_acm_2D_data.m");
    out_file << "positions = " << domain.positions() << ";" << std::endl;
    out_file << "u = " << u_2 << ";" << std::endl;
    out_file.close();
    prn(domain.size());
    double Re = rho / mu;
    prn(Re);
    // ---------
    return 0;
}
