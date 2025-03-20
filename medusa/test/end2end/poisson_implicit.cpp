#include <medusa/Medusa_fwd.hpp>
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"

#include "gtest/gtest.h"

namespace mm {

class PoissonImplicitTest : public ::testing::Test {
  public:
    // Parameters

    Vec2d size;
    Vec2d shift;

    int N;

  protected:
    void SetUp() override {
        size = {1, 1};
        shift = {0, 0};
    }

    template <class BasisType, class WeightType>
    double solve(const BasisType& basis, int support_size,
                 int num_nodes_x, const WeightType& weight) {
        // Closed form solution -- refer to cantilever_beam.nb for reference.
        std::function<double(Vec2d)> analytical = [=](const Vec2d& p) {
            return std::sin(PI*(p[0]-shift[0])/size[0])*std::sin(PI*(p[1]-shift[1])/size[1]);
        };

        // Domain definition
        Vec2d low = shift, high = shift+size;
        BoxShape<Vec2d> box(low, high);
        DomainDiscretization<Vec2d> domain = box.discretizeWithStep(size[0] / num_nodes_x);
        N = domain.size();

        domain.findSupport(FindClosest(support_size));

        WLS<BasisType, WeightType, ScaleToFarthest> approx(basis, weight);
        auto storage = domain.computeShapes<sh::lap>(approx, domain.interior());

        Eigen::SparseMatrix<double, Eigen::RowMajor> M(N, N);
        M.reserve(storage.supportSizes());

        Eigen::VectorXd rhs(N); rhs.setZero();
        auto op = storage.implicitOperators(M, rhs);
        // Set equation on interior
        for (int i : domain.interior()) {
            auto p = domain.pos(i);
            op.lap(i) = -PI*PI*(size.squaredNorm())/size[0]/size[0]/size[1]/size[1]*
                    std::sin(PI*(p[0]-shift[0])/size[0])*std::sin(PI*(p[1]-shift[1])/size[1]);
        }
        for (int i : domain.boundary()) {
            op.value(i) = 0;
        }
        M.makeCompressed();

        Eigen::SparseMatrix<double> M2 = M;
        Eigen::SparseLU<decltype(M2)> solver;
        solver.compute(M2);
        Eigen::VectorXd sol = solver.solve(rhs);

        Range<double> error(N);
        double maxan = 0;
        for (int i = 0; i < N; ++i) {
            double an = analytical(domain.pos(i));
            if (an > maxan) maxan = an;
            error[i] = std::abs(an - sol[i]);
        }
        double L_inf_error = *std::max_element(error.begin(), error.end());
        return L_inf_error / maxan;
    }
};

TEST_F(PoissonImplicitTest, Monomials5NoWeight) {
    Monomials<Vec2d> mon5({{0, 0}, {0, 1}, {1, 0}, {2, 0}, {0, 2}});
    int ss = 5;
    EXPECT_LT(solve(mon5, ss, 29, NoWeight<Vec2d>()), 1e-3);
    EXPECT_LT(solve(mon5, ss, 93, NoWeight<Vec2d>()), 1e-4);
//    EXPECT_LT(solve(mon5, ss, 300, NoWeight<Vec2d>()), 1e-5);
}

TEST_F(PoissonImplicitTest, Monomials5NoWeightShifted) {
    shift = {-123.234, 234.98713};
    Monomials<Vec2d> mon5({{0, 0}, {0, 1}, {1, 0}, {2, 0}, {0, 2}});
    int ss = 5;
    EXPECT_LT(solve(mon5, ss, 29, NoWeight<Vec2d>()), 1e-3);
}
TEST_F(PoissonImplicitTest, Monomials5NoWeightScaledUp) {
    size = 1000.2345;
    Monomials<Vec2d> mon5({{0, 0}, {0, 1}, {1, 0}, {2, 0}, {0, 2}});
    int ss = 5;
    EXPECT_LT(solve(mon5, ss, 29, NoWeight<Vec2d>()), 1e-3);
}
TEST_F(PoissonImplicitTest, Monomials5NoWeightScaledDown) {
    size = 0.00112;
    Monomials<Vec2d> mon5({{0, 0}, {0, 1}, {1, 0}, {2, 0}, {0, 2}});
    int ss = 5;
    EXPECT_LT(solve(mon5, ss, 29, NoWeight<Vec2d>()), 1e-3);
}

TEST_F(PoissonImplicitTest, Gaussians5NoWeight) {
    double sigmaB = 100;
    Gaussians<Vec2d> g5(5, sigmaB);
    int ss = 5;
    EXPECT_LT(solve(g5, ss, 29, NoWeight<Vec2d>()), 1e-3);
    EXPECT_LT(solve(g5, ss, 93, NoWeight<Vec2d>()), 1e-5);  // stagnation error reached
    double err = solve(g5, ss, 167, NoWeight<Vec2d>());
    EXPECT_LT(err, 1e-4);  // error starts increasing
    EXPECT_GT(err, 1e-5);
}

TEST_F(PoissonImplicitTest, Gaussians5Weight) {
    double sigmaB = 100;
    int ss = 5;
    Gaussians<Vec2d> g5(ss, sigmaB);
    GaussianWeight<Vec2d> w(1000.0);
    EXPECT_LT(solve(g5, ss, 29, NoWeight<Vec2d>()), 1e-3);
}

TEST_F(PoissonImplicitTest, Monomials9NoWeight) {
    Monomials<Vec2d> mon9 = Monomials<Vec2d>::tensorBasis(2);
    int ss = 9;
    EXPECT_LT(solve(mon9, ss, 29, NoWeight<Vec2d>()), 1e-3);
    EXPECT_LT(solve(mon9, ss, 93, NoWeight<Vec2d>()), 1e-4);
//    EXPECT_LT(solve(mon9, ss, 300, NoWeight<Vec2d>()), 1e-5);
}

TEST_F(PoissonImplicitTest, Monomials9Weight) {
    Monomials<Vec2d> mon9 = Monomials<Vec2d>::tensorBasis(2);
    GaussianWeight<Vec2d> w(1000.0);  // effectively constant w
    int ss = 9;
    EXPECT_LT(solve(mon9, ss, 29, w), 1e-3);
    EXPECT_LT(solve(mon9, ss, 93, w), 1e-4);
}

TEST_F(PoissonImplicitTest, Gaussians9NoWeight) {
    double sigmaB = 100;
    Gaussians<Vec2d> g9(9, sigmaB);
    int ss = 9;
    EXPECT_LT(solve(g9, ss, 22, NoWeight<Vec2d>()), 1e-2);
    EXPECT_LT(solve(g9, ss, 66, NoWeight<Vec2d>()), 1e-3);
    EXPECT_LT(solve(g9, ss, 149, NoWeight<Vec2d>()), 1e-4);
}

TEST_F(PoissonImplicitTest, Monomials6Weight) {
    GaussianWeight<Vec2d> w(1.0);
    Monomials<Vec2d> mon6(2);
    int ss = 9;
    EXPECT_LT(solve(mon6, ss, 18, w), 1e-2);
    EXPECT_LT(solve(mon6, ss, 59, w), 1e-3);
//    EXPECT_LT(solve(mon6, ss, 188, w), 1e-4);
}

// TEST_F(PoissonImplicitTest, GaussiansWeight) {
//     std::vector<int> ns = {11, 12, 13, 15, 16, 18, 21, 23, 26, 29, 33, 37, 41, 46, 52, 59, 66,
//                            74, 83, 93, 105, 118, 133, 149, 167, 188, 212, 238, 267, 300};
//     std::vector<double> errs;
//     std::vector<int> Ns;
//     for (int n : ns) {
//         double err = solve(mon, ss, n, w);
//         prn(err);
//         errs.push_back(err);
//         Ns.push_back(N);
//     }
//
//     HDF hdf("tmp.h5", HDF::DESTROY);
//     hdf.writeDoubleArray("errs", errs);
//     hdf.writeDoubleArray("Ns", Ns);
// }

}  // namespace mm
