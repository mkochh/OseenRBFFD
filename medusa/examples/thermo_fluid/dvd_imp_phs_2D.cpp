#include <medusa/Medusa_fwd.hpp>
#include <medusa/bits/domains/GeneralFill.hpp>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>

#define EPS 0.0000001  // small number

/// Basic medusa example -- Natural convection on irregular domain - extension of de Vahl Davis test
/// -- http://e6.ijs.si/medusa/wiki/index.php/De_Vahl_Davis_natural_convection_test
/// - implicit RBF-FD solution with PHS basis on scattered nodes
/// - explicit pressure correction http://e6.ijs.si/medusa/wiki/index.php/Fluid_Mechanics#Explicit_Pressure_correction
/// - irregular domain
/// - no internal iterations
/// - using ghost nodes http://e6.ijs.si/medusa/wiki/index.php/Ghost_nodes_(theory)

using namespace mm;  // NOLINT
typedef Vec2d vec_t;
typedef double scal_t;

int main() {
    // case setup
    scal_t dt = 0.01;  // time step
    scal_t time = 2;  // time
    int n = 25;  // no. of support nodes
    scal_t dl = 0.02;  // target distance between nodes
    scal_t mu = 0.071;  // viscosity
    scal_t rho = 1.0;  // density
    scal_t lam = 0.01;  // thermal conductivity
    scal_t c_p = 1.0;  // specificheat capacity
    scal_t g_0 = -500;  // gravitational acceleration
    Vec2d g{0, g_0};
    scal_t beta = 0.71;  // thermal expansion coefficient
    scal_t T_ref = 0;  // reference temperature
    scal_t T_cold = -1.0;  // boundary conditions
    scal_t T_hot = 1.0;

    double Ra = std::abs(g_0) * beta * rho * rho * c_p * std::pow(1, 3) *
                std::abs(T_hot - T_cold) / (lam * mu);
    std::cout << "Pr = " << mu * c_p / lam << " Ra = " << Ra << std::endl;

    //  domain
    BoxShape<vec_t> box(0.0, 1.0);;
    DomainDiscretization<vec_t> domain = box.discretizeBoundaryWithStep(dl);
    //  Obstacles
    Range<double> o_r = {0.24, 0.06, 0.15, 0.05, 0.2, 0.09, 0.09, 0.13};
    Range<vec_t> o_c = {{0.95, 0.13}, {0.02, 0.5}, {0.55, 0.95}, {0.73, 0.4},
                        {0.03, 0.0}, {0.39, 0.62}, {1, 1}, {1, 0.4}};
    Range<int> o_t =    {-16, -15, -14, -15, -14, -16, -14, -16};
    //  cold side obstacles - 16, hot side obstacles - 15, Neumann side -14
    for (auto i=0; i < o_t.size() ; ++i) {
        BallShape<vec_t> ball(o_c[i], o_r[i]);
        DomainDiscretization<vec_t> obstacle = ball.discretizeBoundaryWithStep(dl);
        obstacle.types() = o_t[i];
        domain.subtract(obstacle);
    }
    //  Fill domain
    GeneralFill<vec_t> fill;
    fill(domain, dl);
    //  remove corner nodes
    Range<int> corner = domain.positions().filter([](const vec_t& p) {
        return ((p[0] < EPS && (p[1] < EPS || p[1] > 1 - EPS)) ||
                (p[0] > 1 - EPS && (p[1] < EPS || p[1] > 1 - EPS)));
    });
    domain.removeNodes(corner);
    //  boundary maps
    Range<int> hot_side_obstacles = domain.types() == -15;
    Range<int> cold_side_obstacles = domain.types() == -16;
    Range<int> neumann_side_obstacles = domain.types() == -14;
    Range<int> left_idx = domain.types() == -1;
    Range<int> right_idx = domain.types() == -2;
    Range<int> top_idx = domain.types() == -3;
    Range<int> bottom_idx = domain.types() == -4;
    top_idx = top_idx + neumann_side_obstacles;
    Range<int> interior = domain.interior();
    Range<int> boundary = domain.boundary();
    //  ghost nodes
    Range<int> gh(domain.size(), -100);
    for (int i : boundary) {
        gh[i] = domain.addInternalNode(domain.pos(i) + domain.normal(i) * dl, 9);
    }
    int N = domain.size();
    prn(domain.size());
    domain.findSupport(FindClosest(n));

    //  operator discretization
    Polyharmonic<double, 3> ph;  // construct Polyharmonic
    RBFFD<Polyharmonic<double, 3>, vec_t, ScaleToClosest> appr(ph, Monomials<vec_t>(2));

    //  field initialization
    VectorField2d u_1(domain.size());
    ScalarFieldd p = ScalarFieldd::Zero(domain.size());
    ScalarFieldd T_1(domain.size());
    u_1.setZero();
    T_1.setZero();
    //  dirichlet boundary conditions
    T_1 = (T_cold + T_hot) / 2;
    T_1[left_idx] = T_cold;
    T_1[right_idx] = T_hot;

    ScalarFieldd T_2 = T_1;
    VectorField2d u_2 = u_1;

    //  pressure correction matrix -- note N + 1, for additional constraint
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver_p;
    Eigen::SparseMatrix<double, Eigen::RowMajor> M(N + 1, N + 1);
    Eigen::VectorXd rhs(N + 1);

    auto storage = domain.computeShapes<sh::lap | sh::grad>(appr);
    auto op_v = storage.explicitVectorOperators();  // vector operators
    auto op_s = storage.explicitOperators();  // scalar operators
    auto op_i = storage.implicitOperators(M, rhs);
    Range<int> per_row(N+1, n+1);
    per_row[N] = N;
    M.reserve(per_row);

    //  pressure correction
    for (int i : interior) {
        op_i.lap(i).eval(1);
    }
    for (int i : boundary) {
        op_i.neumann(i, domain.normal(i)).eval(1);
    }
    //  ghost nodes
    for (int i : boundary) {
        op_i.lap(i, gh[i]).eval(1);
    }

    //  regularization
    //  set the last row and column of the matrix
    for (int i = 0; i < N; ++i) {
        M.coeffRef(N, i) = 1;
        M.coeffRef(i, N) = 1;
    }
    //  set the sum of all values
    rhs[N] = 0.0;
    M.makeCompressed();
    Eigen::SparseMatrix<double> MM(M);
    solver_p.compute(MM);
    //  time loop
    for (int time_step = 0; time_step < std::floor(time / dt); ++time_step) {
        prn(time_step / time * dt);
        //  implicit navier stokes
        Eigen::SparseMatrix<double, Eigen::RowMajor> M_v(2 * N, 2 * N);
        Eigen::VectorXd rhs_vec(2 * N);
        Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor>,
                Eigen::IncompleteLUT<double>> solver_u;
        rhs_vec.setZero();
        Range<int> per_row_v(2 * N, n);
        auto op_iv = storage.implicitVectorOperators(M_v, rhs_vec);
        M_v.reserve(per_row_v);
        for (int i : domain.all()) {
            1 / dt * op_iv.value(i) + op_iv.grad(i, u_1[i]) + (-1 * mu / rho) * op_iv.lap(i)
                    = u_1[i] / dt + 1 * g * (1 - beta * (T_1[i] - T_ref));
        }
        M_v.makeCompressed();
        solver_u.compute(M_v);
        Eigen::VectorXd solution = solver_u.solveWithGuess(rhs_vec, u_1.asLinear());
        u_2 = VectorField2d::fromLinear(solution);
        //  pressure correction
        for (int i : interior) rhs(i) = rho / dt * op_v.div(u_2, i);
        for (int i : boundary) rhs(i) = rho / dt * u_2[i].dot(domain.normal(i));
        //  ghost nodes
        for (int i : boundary) {
            rhs(gh[i]) = 0;
        }
        ScalarFieldd P_c = solver_p.solve(rhs).head(N);
        for (int i = 0; i < interior.size(); ++i) {
            int c = interior[i];
            u_2[c] -= dt / rho * op_s.grad(P_c, c);
        }
        //  force BCs
        u_2[boundary] = 0;
        //  ghost nodes
        for (int i : boundary) {
            u_2(gh[i]).setZero();
        }
        //  implicit heat transfer
        Eigen::SparseMatrix<double, Eigen::RowMajor> M_T(N, N);
        Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor>,
                Eigen::IncompleteLUT<double>> solver_T;
        Eigen::VectorXd rhs_T(N);
        rhs_T.setZero();
        auto op_t = storage.implicitOperators(M_T, rhs_T);

        Range<int> per_row(N, n);
        M_T.reserve(per_row);
        //  heat transfer
        for (int i : interior) {
            op_t.value(i) + (-dt * lam / rho / c_p) * op_t.lap(i) +
            dt * op_t.grad(i, u_2[i]) = T_1(i);
        }
        //  BCs
        for (int i : top_idx + bottom_idx) {
            op_t.neumann(i, domain.normal(i)) = 0;
        }
        for (int i : left_idx) op_t.value(i) = T_hot;
        for (int i : right_idx) op_t.value(i) = T_cold;
        for (int i : cold_side_obstacles) op_t.value(i) = T_cold;
        for (int i : hot_side_obstacles) op_t.value(i) = T_hot;
        //  ghost nodes
        for (int i : top_idx + bottom_idx) {
            op_t.value(i, gh[i]) + (-dt * lam / rho / c_p) * op_t.lap(i, gh[i]) +
            dt * op_t.grad(i, u_2[i], gh[i]) = T_1(i);
        }
        for (int i : left_idx) op_t.value(gh[i]) = T_hot;
        for (int i : right_idx) op_t.value(gh[i]) = T_cold;
        for (int i : hot_side_obstacles) op_t.value(gh[i]) = T_hot;
        for (int i : cold_side_obstacles) op_t.value(gh[i]) = T_cold;

        M_T.makeCompressed();
        solver_T.compute(M_T);
        T_2 = solver_T.solveWithGuess(rhs_T, T_1);

        T_1.swap(T_2);
        u_1.swap(u_2);
    }

    std::ofstream out_file("dvd_imp_phs_2D_data.m");
    out_file << "positions = " << domain.positions() << ";" << std::endl;
    out_file << "u = " << u_2 << ";" << std::endl;
    out_file << "T = " << T_2 << ";" << std::endl;
    out_file.close();

    return 0;
}

