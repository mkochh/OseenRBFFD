//
// Created by blaz on 5/20/20.
//
//
// Benchmark case for FDFD in 2D, scattering from a dielectric cylinder for a TE planewave.
// The analytical solution is reasonably easy to find in terms of Hankel functions.
// The solution uses the radial PML variant

#include <medusa/Medusa_fwd.hpp>
#include <medusa/bits/domains/BasicRelax.hpp>
#include <medusa/bits/domains/GeneralFill.hpp>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>
#include <algorithm>
#include <complex>

using namespace mm; // NOLINT

constexpr std::complex<double> operator"" _i(long double d) {
    return std::complex<double>{0.0, static_cast<double>(d)};
}

int main() {
    XML conf("dielectric_cylinder.xml");
    std::string outfile = conf.get<std::string>("IO.outfile");

    // solver
    int ffact = conf.get<int>("solver.ffact");
    int maxiter = conf.get<int>("solver.maxiter");
    double dtol = conf.get<double>("solver.droptol");
    double soltol = conf.get<double>("solver.errtol");

    // Physics
    double mu0 = 1.25663706212*1e-6;  // [H/m = N/A^2]
    double eps0 = 8.8541878128*1e-12;  // [F/m]
    double c0 = 1.0/sqrt(mu0*eps0);
    double epsr = conf.get<double>("case.epsr");
    double wavelength = conf.get<double>("fourier.wavelength");
    double k0 = 2*PI/wavelength;
    double omega = k0*c0;

    // Geometry
    double r = conf.get<double>("geometry.r");
    double a = conf.get<double>("geometry.a");
    double tol = 1e-5;
    double t0 = conf.get<double>("case.theta0");

    std::cout << "Solving for:" <<
              "\nk0 = " << k0 <<
              "\nomega = " << omega <<
              "\nwavelength = " << wavelength << std::endl;

    // Approximation
    int supp = conf.get<int>("approx.supp");
    int num_monoms = conf.get<int>("approx.num_monoms");

    // Fill
    double dx_cyl = conf.get<double>("case.dx_cylinder");
    double dx_omega = conf.get<double>("case.dx_omega");
    double alpha = conf.get<double>("case.alpha");
    double power = conf.get<double>("case.power");
    auto dx = [=](const Vec2d &p){return (dx_cyl-dx_omega)/
                          (1.0+std::pow(std::sqrt(p[0]*p[0]+p[1]*p[1])/alpha, power))+dx_omega;};

    // PML
    double d = conf.get<double>("PML.d");
    double m = conf.get<double>("PML.m");
    double sigma_max = conf.get<double>("PML.sigma");

    BallShape <Vec2d> outbox(0, r);  // outer
    BallShape <Vec2d> cylinder(0, a);  // inner
    auto outer_shape = outbox - cylinder;

    // fill outer boundary
    DomainDiscretization<Vec2d> cyl_dom = cylinder.discretizeBoundaryWithDensity(dx);
    DomainDiscretization<Vec2d> outer_dom = (outbox-cylinder).discretizeBoundaryWithDensity(dx);

    // Node bookkeeping
    // delete nodes of outer domain on shared boundary
    Range<int> to_del = outer_dom.positions().filter([=]
            (const Vec2d& p){ return ((p).norm()-a < tol);});
    outer_dom.removeNodes(to_del);

    // make index arrays
    Range<int> interface_c_idx(0);  // cylinder index
    Range<int> interface_o_idx(0);  // outer index

    for (int i : cyl_dom.boundary()) {
        int id = outer_dom.addBoundaryNode(cyl_dom.pos(i), -6, -1 * cyl_dom.normal(i));
        interface_c_idx.push_back(i);
        interface_o_idx.push_back(id);
    }

    // outer boundary
    Range<int> outer_bnd = outer_dom.positions().filter([=]
            (const Vec2d& p) {return r-p.norm() < tol;});

    // fill domains
    GeneralFill<Vec2d> fill;
    fill(cyl_dom, dx);
    fill(outer_dom, dx);

    int Ncyl = cyl_dom.size();
    int Nouter = outer_dom.size();

    std::cout << "Domain size = " << Ncyl + Nouter << "\n";

    BasicRelax relax;
    relax(cyl_dom, dx);
    relax(outer_dom, dx);

    // filter for PML
    auto PML = (outer_dom.positions()).filter([=]
            (const Vec2d& p){ return (std::sqrt(r - d) < std::sqrt(p.norm()) &&
            std::sqrt(p.norm()) < std::sqrt(r-tol));});
    auto outer_inter = (outer_dom.positions()).filter([=]
            (const Vec2d& p) { return (std::sqrt(a + tol) < std::sqrt(p.norm())
            && std::sqrt(p.norm()) <= std::sqrt(r-d));});

    // find support
    cyl_dom.findSupport(FindClosest(supp));
    outer_dom.findSupport(FindClosest(supp));

    // make rbffd engine
    Polyharmonic<double, 3> phs;
    RBFFD<Polyharmonic<double, 3>, Vec2d, ScaleToClosest,
    Eigen::PartialPivLU<Eigen::MatrixXd>> rbffd(phs, num_monoms);

    // make M and rhs
    Eigen::SparseMatrix <std::complex<double>, Eigen::RowMajor> M(Ncyl + Nouter, Ncyl + Nouter);
    Eigen::VectorXcd rhs(Ncyl + Nouter); rhs.setZero();

    // make operators
    auto storage_inner = cyl_dom.computeShapes<sh::lap | sh::d1 | sh::d2>(rbffd);
    auto storage_outer = outer_dom.computeShapes<sh::lap | sh::d1 | sh::d2>(rbffd);

    auto op_inner = storage_inner.implicitOperators(M, rhs);
    auto op_outer = storage_outer.implicitOperators(M, rhs);

    // offset operator write
    op_inner.setRowOffset(Nouter);
    op_inner.setColOffset(Nouter);

    // reserve storage
    auto store = storage_inner.supportSizes() + storage_outer.supportSizes();
    Range<int> reserve = store;
    for (int& i : reserve) i *= 2;
    M.reserve(reserve);

    ExplicitOperators<decltype(storage_outer)> exop(storage_outer);
    Eigen::VectorXcd sw(Nouter);

    for (int i : outer_inter) {
        sw[i] = 1.0;
    }
    for (int i : PML) {
        Vec2d pos = outer_dom.pos(i);
        double rad = std::sqrt(std::pow(pos[0], 2) + std::pow(pos[1], 2));
        sw[i] = 1.0-1.0_i*sigma_max/omega/eps0*std::pow((rad+d-r)/d, m);
    }

    // Construct system
    double l_in = epsr*std::pow(omega/c0, 2.0);
    double l_out = std::pow(omega/c0, 2.0);
    for (int i : cyl_dom.interior()) {
        op_inner.lap(i) + l_in*op_inner.value(i) = 0.0;
    }
    for (int i : outer_inter) {
        op_outer.lap(i) + l_out*op_outer.value(i) = 0.0;
    }
    for (int c = 0; c < interface_c_idx.size(); ++c) {
        int i = interface_c_idx[c];  // cyl index
        int j = interface_o_idx[c];  // outer index
        Vec2d pos = cyl_dom.pos(i);
        double x = pos[0];
        double y = pos[1];

        // calculate normal derivative of the source funciton at the i-th node
        // get normal, must point out of the cylinder
        Vec2d normal = cyl_dom.normal(i);
        std::complex<double> incident = std::exp(1.0_i*k0*(x*std::cos(t0)+y*std::sin(t0)));
        std::complex<double> din_dx =
                1.0_i*k0*std::cos(t0)*std::exp(1.0_i*k0*(x*std::cos(t0)+y*std::sin(t0)));
        std::complex<double> din_dy =
                1.0_i*k0*std::sin(t0)*std::exp(1.0_i*k0*(x*std::cos(t0)+y*std::sin(t0)));
        std::complex<double> dui_dn = normal[0]*din_dx + normal[1]*din_dy;

        // continuity of the fields, and the incident field
        op_inner.value(i) + (-1)*op_outer.value(j, Nouter + i) = incident;

        // continuity of derivatives
        op_outer.neumann(j, outer_dom.normal(j)) +
                (1/epsr)*op_inner.neumann(i, cyl_dom.normal(i), j - Nouter)
                    = dui_dn;
    }
    // SC-PML - where either x or y is constant
    for (int i : PML) {
        1.0/(sw[i]*sw[i])*op_outer.lap(i) + l_out*op_outer.value(i) +
        ((-1.0)/(sw[i]*sw[i]*sw[i])*exop.d1(sw, 0, i)*op_outer.der1(i, 0)
         +(-1.0)/(sw[i]*sw[i]*sw[i])*exop.d1(sw, 1, i)*op_outer.der1(i, 1)) = 0.0;
    }
    for (int i : outer_bnd) {
        op_outer.value(i) = 0.0;
    }

    HDF out(outfile+".h5", HDF::DESTROY);
    out.writeXML("config", conf);
    out.writeDomain("cyl_dom", cyl_dom);
    out.writeDomain("outer_dom", outer_dom);

    Eigen::BiCGSTAB<decltype(M), Eigen::IncompleteLUT<std::complex<double>>> solver;

    // customize solver
    solver.preconditioner().setFillfactor(ffact);
    solver.preconditioner().setDroptol(dtol);
    solver.setMaxIterations(maxiter);
    solver.setTolerance(soltol);
    solver.compute(M);
    Eigen::VectorXcd sol = solver.solve(rhs);
    std::cout << "Iterations: " << solver.iterations() << "\n";
    std::cout << "Error: " << solver.error() << "\n";

    ScalarField<double> rsol = sol.real();
    ScalarField<double> csol = sol.imag();

    // out types for debugging
    out.writeIntArray("pml_nodes", PML);
    out.writeIntArray("outer_inter_nodes", outer_inter);
    out.writeDoubleArray("rsol", rsol);
    out.writeDoubleArray("csol", csol);
    out.writeDoubleAttribute("wavelength", wavelength);
    out.close();

    return 0;
}
