#include <medusa/Medusa_fwd.hpp>
#include <medusa/bits/domains/GeneralFill.hpp>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>


/// Classical 3D Point contact problem.
/// http://e6.ijs.si/medusa/wiki/index.php/Linear_elasticity#Point_contact_3D

using namespace mm;  // NOLINT

int main() {
    XML conf("point_contact3d.xml");

    // Physical parameters
    const double E = conf.get<double>("phy.E");
    const double nu = conf.get<double>("phy.nu");
    const double P = -conf.get<double>("phy.P");

    // Derived parameters
    double lam = E * nu / (1 - 2 * nu) / (1 + nu);
    const double mu = E / 2 / (1 + nu);

    // Closed form solution
    std::function<Vec3d(Vec3d)> analytical = [=](const Vec3d& p) {
        double x = p[0], y = p[1], z = p[2];
        double r = std::sqrt(x*x+y*y);
        double c = x/r, s = y/r;
        double R = p.norm();
        double u = P*r/4/PI/mu * (z/R/R/R - (1-2*nu) / (R*(R+z)));
        double w = P/4/PI/mu * (z*z/R/R/R + 2*(1-nu)/R);
        return Vec3d(u*c, u*s, w);
    };

    // Domain definition
    double a = conf.get<double>("num.a");
    double b = conf.get<double>("num.b");
    Vec3d low = -a, high = -b;
    BoxShape<Vec3d> box(low, high);

    double step = conf.get<double>("num.h");

    std::function<double(Vec3d)> density = [=] (const Vec3d& v) {
        double r = (v-low).head<2>().norm() / (high-low).head<2>().norm();
        return step*(9*(1-r)+1);
    };
    GeneralFill<Vec3d> fill; fill.seed(0);
    DomainDiscretization<Vec3d> domain = box.discretizeWithDensity(density, fill);
    int N = domain.size();
    prn(N);
    domain.findSupport(FindClosest(15));  // the support for each node is the closest 9 nodes

    // Construct the approximation engine, in this case a weighted least squares using monomials as
    // basis functions, no weight, and scale to farthest.
    WLS<Gaussians<Vec3d>, GaussianWeight<Vec3d>, ScaleToFarthest> wls({15, 100.0}, 1.0);

    // Compute the shapes using the defined approximation engine
    auto shapes = domain.computeShapes(wls);
    std::cerr << "Shapes computed, solving system." << std::endl;

    Eigen::SparseMatrix<double, Eigen::RowMajor> M(3*N, 3*N);
    Eigen::VectorXd rhs(3*N); rhs.setZero();
    M.reserve(shapes.supportSizesVec());

    // Construct implicit operators over our storage.
    auto op = shapes.implicitVectorOperators(M, rhs);

    // Set the governing equations and the boundary conditions.
    for (int i : domain.interior()) {
        (lam+mu)*op.graddiv(i) + mu*op.lap(i) = 0.0;
    }
    for (int i : domain.boundary()) {
        op.value(i) = analytical(domain.pos(i));
    }

    Eigen::BiCGSTAB<decltype(M), Eigen::IncompleteLUT<double>> solver;
    solver.preconditioner().setFillfactor(20);
    solver.preconditioner().setDroptol(1e-4);
    solver.setMaxIterations(100);
    solver.setTolerance(1e-14);
    solver.compute(M);
    Eigen::VectorXd solution = solver.solve(rhs);
    prn(solver.iterations());
    prn(solver.error());
    VectorField3d u = VectorField3d::fromLinear(solution);

    auto eop = shapes.explicitVectorOperators();
    VectorField<double, 6> stress(N);
    for (int i = 0; i < N; ++i) {
        Eigen::Matrix3d grad = eop.grad(u, i);
        Eigen::Matrix3d eps = 0.5*(grad + grad.transpose());
        Eigen::Matrix3d s = lam * eps.trace() * Eigen::Matrix3d::Identity(3, 3) + 2*mu*eps;
        stress[i][0] = s(0, 0);
        stress[i][1] = s(1, 1);
        stress[i][2] = s(2, 2);
        stress[i][3] = s(0, 1);
        stress[i][4] = s(0, 2);
        stress[i][5] = s(1, 2);
    }


    // Write the solution into a file.
    std::ofstream out_file("point_contact3d_data.m");
    out_file << "E = " << E << ";\n";
    out_file << "nu = " << nu << ";\n";
    out_file << "P = " << P << ";\n";
    out_file << "positions = " << domain.positions() << ";\n";
    out_file << "displacement = " << u << ";" << ";\n";
    out_file << "stress = " << stress << ";" << std::endl;
    out_file.close();
    return 0;
}
