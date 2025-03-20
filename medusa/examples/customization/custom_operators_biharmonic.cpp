#include <medusa/Medusa.hpp>
#include <medusa/bits/domains/GeneralFill.hpp>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>

using namespace mm;  // NOLINT
using namespace std;  // NOLINT
using namespace Eigen;  // NOLINT

/// Solving Biharmonic equation to demonstrate custom operators.
/// http://e6.ijs.si/medusa/wiki/index.php/Customization

constexpr int dim = 2;  // Change to 1 or 3.

/// Represents the Biharmonic operator
template <int dim>
struct Biharmonic : public Operator<Biharmonic<dim>> {
    typedef Vec<double, dim> vec;
    static std::string type_name() { return format("Biharmonic<%d>", dim); }
    static std::string name() { return type_name(); }

    template <int k>
    double applyAt0(const RBFBasis<Polyharmonic<double, k>, vec>&, int index,
                    const std::vector<vec>& support, double scale) const {
        static_assert(k != -1, "If dynamic k is desired it can be obtained from basis.rbf().");
        double r = support[index].norm();
        return k*(k-2)*(dim+k-2)*(dim+k-4)*ipow<k-4>(r) / ipow<4>(scale);
    }

    double applyAt0(const Monomials<vec>& mon, int idx, const std::vector<vec>& q, double s) const {
        double result = 0;
        std::array<int, dim> orders;
        for (int d = 0; d < dim; ++d) { orders[d] = 0; }
        for (int d = 0; d < dim; ++d) {
            orders[d] = 4;
            result += mon.evalOpAt0(idx, Derivative<dim>(orders), q);
            orders[d] = 0;

            for (int d2 = 0; d2 < d; ++d2) {
                orders[d] = 2;
                orders[d2] = 2;
                result += 2*mon.evalOpAt0(idx, Derivative<dim>(orders), q);
                orders[d] = 0;
                orders[d2] = 0;
            }
        }
        return result / ipow<4>(s);
    }
};

int main() {
    typedef Vec<double, dim> vec;
    BallShape<vec> b(0.0, 1.0);
    double dx = 0.05;
    DomainDiscretization<vec> domain = b.discretizeBoundaryWithStep(dx);
    auto fn = [=](const vec&) { return dx; };
    GeneralFill<vec> fill; fill.seed(1337);
    fill(domain, fn);

    auto gh = domain.addGhostNodes(dx);

    int N = domain.size();
    Monomials<vec> mon(2);
    int n = 2*mon.size();
    domain.findSupport(FindClosest(n));
    RBFFD<Polyharmonic<double, 5>, vec, ScaleToClosest> approx({}, mon);

    std::tuple<Biharmonic<dim>, Der1s<dim>> ops;
    cout << ops << endl;
    RaggedShapeStorage<vec, decltype(ops)> storage;
    storage.resize(domain.supportSizes());
    computeShapes(domain, approx, domain.interior(), std::tuple<Biharmonic<dim>>(), &storage);
    computeShapes(domain, approx, domain.boundary(), std::tuple<Der1s<dim>>(), &storage);

    Eigen::SparseMatrix<double, Eigen::RowMajor> M(N, N);
    M.reserve(Range<int>(N, n));
    Eigen::VectorXd rhs(N); rhs.setZero();
    auto op = storage.implicitOperators(M, rhs);
    for (int i : domain.interior()) {
        op.apply<Biharmonic<dim>>(i) = 0.0;
    }
    for (int i : domain.boundary()) {
        op.neumann(i, domain.normal(i)) = 1.0;
        op.value(i, gh[i]) = 1.0;
    }

    SparseLU<decltype(M)> solver;
    solver.compute(M);
    ScalarFieldd u = solver.solve(rhs);

    std::ofstream out_file("custom_operators_biharmomic_data.m");
    out_file << "positions = " << domain.positions() << ";" << std::endl;
    out_file << "solution = " << u << ";" << std::endl;
    out_file.close();

    return 0;
}

