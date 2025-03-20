#include <medusa/Medusa.hpp>
#include <medusa/bits/domains/GeneralFill.hpp>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <iostream>
#include <sstream>
#include "gtest/gtest.h"

namespace mm {
constexpr int dim = 2;  // Change to 1 or 3.


TEST(Approximations, printOperators) {
    Der1s<3> der1;
    Der2s<3> der2;
    Lap<3> lap;
    std::stringstream out;
    out << der1;
    EXPECT_EQ(out.str(), "Der1s<3>");
    out.str("");
    out.clear();
    out << der2;
    EXPECT_EQ(out.str(), "Der2s<3>");
    out.str("");
    out.clear();
    out << lap;
    EXPECT_EQ(out.str(), "Laplacian<3>");
}

/// [Biharmonic def]
/// Represents the Biharmonic operator
template <int dim>
struct Biharmonic : public Operator<Biharmonic<dim>> {
    typedef Vec<double, dim> vec;
    static std::string type_name() { return format("Biharmonic<%d>", dim); }

    /// Implement only the functions you need.
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
/// [Biharmonic def]


template <int dim>
int find_idx(const Eigen::Matrix<int, dim, Eigen::Dynamic>& powers,
             const Eigen::Matrix<int, dim, 1>& mon) {
    int idx = -1;
    for (int i = 0; i < powers.cols(); ++i) {
        EXPECT_FALSE(powers.col(i) == mon && idx != -1);
        if (powers.col(i) == mon) {
            idx = i;
        }
    }
    EXPECT_NE(idx, -1);
    return idx;
}

TEST(Approximations, BiharmonicCostumOp) {
    /// [Biharmonic usage]
    typedef Vec<double, dim> vec;
    BallShape<vec> b(0.0, 1.0);
    double dx = 0.05;
    DomainDiscretization<vec> domain = b.discretizeBoundaryWithStep(dx);
    auto fn = [=](const vec&) { return dx; };
    GeneralFill<vec> fill; fill.seed(1337);
    fill(domain, fn);

    Monomials<vec> mon(6);
    Polyharmonic<double, 5> rbf;
    int n = 2*mon.size();
    domain.findSupport(FindClosest(n));
    RBFFD<Polyharmonic<double, 5>, vec, ScaleToClosest> approx(rbf, mon);
    RBFBasis<Polyharmonic<double, 5>, Vec2d> basis(n);

    Biharmonic<2> bih;
    double val1 = mon.evalOpAt0(0, bih);
    double val2 = mon.evalOpAt0(mon.size()-1, bih);
    double val3 = basis.evalOpAt0(0, bih, domain.supportNodes(0));

    approx.compute(domain.pos(0), domain.supportNodes(0));
    auto shape = approx.getShape(bih);
    /// [Biharmonic usage]

    EXPECT_EQ(0, val1);
    EXPECT_EQ(0, val2);
    EXPECT_EQ(225, val3);
    (void) shape;

    int x4 = find_idx(mon.powers(), {4, 0});
    int x2y2 = find_idx(mon.powers(), {2, 2});
    EXPECT_EQ(24, mon.evalOpAt0(x4, bih));
    EXPECT_EQ(8, mon.evalOpAt0(x2y2, bih));
}

}  // namespace mm

