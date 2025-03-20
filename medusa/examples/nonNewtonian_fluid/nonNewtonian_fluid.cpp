#include <medusa/Medusa.hpp>
#include <medusa/bits/domains/BasicRelax.hpp>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>

#define EPS 1e-10  // a small number

using namespace mm;  // NOLINT
constexpr int dim = 2;
typedef double scal_t;
typedef Vec<scal_t, dim> vec_t;

int main() {
    // case setup
    scal_t dt = 1e-3;  // time step
    scal_t spatial_step = 0.02;  // spatial discretization
    scal_t end_time = 2;  // dimensional end time
    scal_t power_index = 0.7;  // non-Newtonian power index

    scal_t mu = 0.071;  // viscosity
    scal_t rho = 1.0;  // density
    scal_t lam = 0.01;  // thermal conductivity
    scal_t c_p = 1.0;  // specific heat capacity
    scal_t g_0 = -500;  // gravitational acceleration
    scal_t beta = 0.71;  // thermal expansion coefficient
    scal_t T_ref = 0;  // reference temperature
    scal_t T_left = -1.0;  // left boundary temperature
    scal_t T_right = 1.0;  // left boundary temperature
    scal_t height = 1.0;  // cavity height

    int n = 12;  // number of support nodes
    int m = 2;  // monomial basis
    scal_t sigma_w = 1;  // weight function weight

    double alpha = lam / c_p / rho;

    double ra = rho * beta * std::abs(g_0 * (T_right - T_left)) * std::pow(height, power_index + 2)
                / std::pow(alpha, power_index) / mu;
    prn(ra)
    double pr = mu * std::pow(alpha, power_index - 2) * std::pow(height, 1 - power_index) / rho;
    prn(pr)

    BoxShape<vec_t> box({0, 0}, {1, height});
    DomainDiscretization<vec_t> domain = box.discretizeBoundaryWithStep(spatial_step);
    domain = box.discretizeWithStep(spatial_step);

    // remove corner nodes (improves stability at the price of a small accuracy loss)
    Range<int> corner = domain.positions().filter([&](const vec_t& p) {
        return ((p[0] < EPS && (p[1] < EPS || p[1] > height - EPS)) ||
                (p[0] > 1 - EPS && (p[1] < EPS || p[1] > height - EPS)));
    });
    domain.removeNodes(corner);

    // boundary maps
    Range<int> top_idx = domain.positions().filter([&](const vec_t& p) {
        return (p[1] >  height - EPS) && ((p[0] > EPS) && (p[0] < 1 - EPS)); });
    Range<int> left_idx = domain.positions().filter([&](const vec_t& p) {
        return p[0] < EPS; });
    Range<int> bottom_idx = domain.positions().filter([&](const vec_t& p) {
        return (p[1] <  EPS) && ((p[0] > EPS) && (p[0] < 1 - EPS)); });
    Range<int> right_idx = domain.positions().filter([&](const vec_t& p) {
        return p[0] > 1 - EPS; });

    Range<int> interior = domain.interior();
    Range<Range<int>> edges = {top_idx, bottom_idx, left_idx, right_idx};

    int N = domain.size();

    WLS<Monomials<vec_t>, GaussianWeight<vec_t>, ScaleToClosest> appr(m, sigma_w);

    // only interior nodes in edge node support to reduce issues with Neumann BCs
    domain.findSupport(FindClosest(n).forNodes(domain.interior()));
    domain.findSupport(FindClosest(n).forNodes(domain.boundary())
                                        .searchAmong(domain.interior())
                                        .forceSelf(true));

    // field declarations
    Vec2d g{0, g_0};
    VectorField2d u_1(domain.size());
    ScalarFieldd p = ScalarFieldd::Zero(domain.size());
    ScalarFieldd T_1(domain.size());
    ScalarFieldd viscosity(domain.size());
    Range<Eigen::Matrix<scal_t, dim, dim>> v_grad(domain.size());
    u_1.setZero();
    T_1.setZero();
    T_1 = T_ref;
    T_1[left_idx] = T_left;
    T_1[right_idx] = T_right;

    ScalarFieldd T_2 = T_1;
    VectorField2d u_2 = u_1;

    // operators
    auto storage = domain.computeShapes<sh::lap | sh::grad>(appr);
    auto op_v = storage.explicitVectorOperators();  // vector operators
    auto op_s = storage.explicitOperators();  // scalar operators

    // pressure correction matrix -- note N + 1, for additional constraint
    Eigen::SparseLU<Eigen::SparseMatrix<scal_t>> solver_p;

    Eigen::SparseMatrix<scal_t, Eigen::RowMajor> M(N + 1, N + 1);
    Eigen::VectorXd rhs(N + 1);
    rhs.setZero();
    // construct implicit operators over our storage
    auto op_i = storage.implicitOperators(M, rhs);
    // pressure correction
    for (int i : domain.interior()) {
        op_i.lap(i).eval(1);
    }
    for (int i : domain.boundary()) {
        op_i.neumann(i, domain.normal(i)).eval(1);
    }
    // regularization, set the last row and column of the matrix
    for (int i = 0; i < N; ++i) {
        M.coeffRef(N, i) = 1;
        M.coeffRef(i, N) = 1;
    }
    // set the sum of all values
    rhs[N] = 0.0;
    M.makeCompressed();
    Eigen::SparseMatrix<scal_t, Eigen::ColMajor> MM(M);
    solver_p.compute(M);
    if (solver_p.info() != Eigen::Success) {
        std::cout << "LU factorization failed with error:" << solver_p.lastErrorMessage()
                  << std::endl;
    }

// Most inner loops can be parallelized with OpenMP

    // time loop
    for (int time_step = 0; time_step < std::floor(end_time / dt); ++time_step) {
        int i;
        for (i = 0; i < domain.size(); ++i) {
            v_grad[i] = op_v.grad(u_1, i);
            double norm = ((v_grad[i] + v_grad[i].transpose())).squaredNorm() / 2;
            // limit minimal norm value to avoid infinite viscosity
            norm = std::max(norm, EPS);
            viscosity[i] = mu * std::pow(norm, (power_index - 1) / 2);
        }

        for (i = 0; i < interior.size(); ++i) {
            int c = interior[i];
            u_2[c] = u_1[c] + dt * ((v_grad[c] * op_s.grad(viscosity, c)
                                     + viscosity[c] * op_v.lap(u_1, c)) / rho
                                    - v_grad[c] * u_1[c]
                                    - g * (beta * (T_1[c] - T_ref)));
        }


        // pressure correction
        for (i = 0; i < interior.size(); ++i) {
            int c = interior[i];
            rhs(c) = rho / dt * op_v.div(u_2, c);
        }
        for (i = 0; i < edges.size(); ++i) {
            for (int j : edges[i]) rhs(j) = rho / dt * u_2[j].dot(domain.normal(j));
        }

        p = solver_p.solve(rhs).head(N);

        for (i = 0; i < interior.size(); ++i) {
            int c = interior[i];
            u_2[c] -= dt / rho * op_s.grad(p, c);
        }

        // explicit heat transfer
        for (i = 0; i < interior.size(); ++i) {
            int c = interior[i];
            T_2[c] = T_1[c] + dt * alpha * op_s.lap(T_1, c)
                     - dt * u_2[c].transpose() * op_s.grad(T_1, c);
        }

        // can cause issues if top and bottom edges do not have the same amount of nodes
        for (i = 0; i < top_idx.size(); ++i) {
            int c = top_idx[i];
            T_2[c] = op_s.neumann(T_2, c, vec_t{0, -1}, 0.0);
            c = bottom_idx[i];
            T_2[c] = op_s.neumann(T_2, c, vec_t{0, 1}, 0.0);
        }

        T_1.swap(T_2);
        u_1.swap(u_2);
    }

    std::ofstream out_file("nonNewtonian_fluid_data.m");
    out_file << "positions = " << domain.positions() << ";" << std::endl;
    out_file << "u = " << u_2 << ";" << std::endl;
    out_file << "T = " << T_2 << ";" << std::endl;
    out_file << "Ra = " << ra << ";" << std::endl;
    out_file << "Pr = " << pr << ";" << std::endl;
    out_file << "power_index = " << power_index << ";" << std::endl;
    out_file.close();

// ---------
    return 0;
}
