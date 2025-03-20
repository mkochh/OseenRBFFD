#include <medusa/bits/approximations/Operators.hpp>
#include <medusa/bits/approximations/RBFBasis.hpp>
#include <medusa/bits/approximations/Gaussian.hpp>
#include <medusa/bits/types/Vec.hpp>
#include <medusa/bits/types/Range.hpp>

#include "gtest/gtest.h"

namespace mm {

TEST(Approximations, RBFBasis3D) {
    double s = 0.4;
    Gaussian<double> g(s);
    RBFBasis<Gaussian<double>, Vec3d> basis(3, g);
    Vec3d p(1.33, 4.55, -0.48);
    Range<Vec3d> support = {{1.3, 4.5, -0.5}, {1.2, 4.3, -0.6}, {1.4, 4.6, -0.49}};
    EXPECT_NEAR(0.97652981169681979, basis.eval(0, p, support), 1e-15);
    EXPECT_NEAR(0.55640991454257032, basis.eval(1, p, support), 1e-15);
    EXPECT_NEAR(0.95420666596918832, basis.eval(2, p, support), 1e-15);

    EXPECT_NEAR(-0.36619867938630742, basis.evalOp(0, p, Der1<3>(0), support), 1e-14);
    EXPECT_NEAR(-0.9041661111316767, basis.evalOp(1, p, Der1<3>(0), support), 1e-14);
    EXPECT_NEAR(0.83493083272303978, basis.evalOp(2, p, Der1<3>(0), support), 1e-14);

    EXPECT_NEAR(0.228874174616442138, basis.evalOp(0, p, Der2<3>(0, 1), support), 1e-14);
    EXPECT_NEAR(1.356249166697515158, basis.evalOp(1, p, Der2<3>(0, 2), support), 1e-14);
    EXPECT_NEAR(-0.074547395778842837, basis.evalOp(2, p, Der2<3>(1, 2), support), 1e-14);
//    EXPECT_NEAR(-0.065228971306487483, basis.evalOp(2, p, {1, 1, 1}, support), 1e-14);

    EXPECT_NEAR(-12.069298141440380, basis.evalOp(0, p, Der2<3>(0), support), 1e-12);
    EXPECT_NEAR(-1.5214333600773410, basis.evalOp(1, p, Der2<3>(1), support), 1e-12);
    EXPECT_NEAR(-11.912673845459090, basis.evalOp(2, p, Der2<3>(2), support), 1e-12);

    EXPECT_NEAR(basis.evalOp(1, p, Lap<3>(), support),
                basis.evalOp(1, p, Der2<3>(0), support) +
                basis.evalOp(1, p, Der2<3>(1), support) +
                basis.evalOp(1, p, Der2<3>(2), support),
                1e-14);

    // Checked independently using code in Mathematica below.
    // s = 0.4`20;
    // p = {1.33`20, 4.55`20, -0.48`20};
    // supp = {{1.3`20, 4.5`20, -0.5`20}, {1.2`20, 4.3`20, -0.6`20}, {1.4`20, 4.6`20, -0.49`20}};
    // basis = Table[With[{q = sn, sc = s}, Function[x, Exp[-(x - q).(x - q)/sc^2]]], {sn, supp}]
    // D[basis[[3]][{x, y, z}], {y, 0}, {z, 0}, {x, 2}] /. {x -> p[[1]], y -> p[[2]], z -> p[[3]]}
}

TEST(Approximations, RBFBasis2D) {
    double s = 1.23;
    Gaussian<double> g(s);
    RBFBasis<Gaussian<double>, Vec2d> basis(2, g);
    Vec2d p(0.6, -0.33);
    Range<Vec2d> support = {{-0.77, 0.64}, {0.49, -2.1}};
    EXPECT_NEAR(0.15528149718493840, basis.eval(0, p, support), 1e-15);
    EXPECT_NEAR(0.12508158425064380, basis.eval(1, p, support), 1e-15);

    EXPECT_NEAR(-0.28122896575235060, basis.evalOp(0, p, Der1<2>(0), support), 1e-14);
    EXPECT_NEAR(-0.29267552927971390, basis.evalOp(1, p, Der1<2>(1), support), 1e-14);

    EXPECT_NEAR(-0.3606214512258313188, basis.evalOp(0, p, Der2<2>(0, 1), support), 1e-14);

    EXPECT_NEAR(0.050053899223872840, basis.evalOp(0, p, Der2<2>(1), support), 1e-12);
    EXPECT_NEAR(-0.16270845136299020, basis.evalOp(1, p, Der2<2>(0), support), 1e-12);

    EXPECT_NEAR(basis.evalOp(1, p, Lap<2>(), support),
                basis.evalOp(1, p, Der2<2>(1), support) + basis.evalOp(1, p, Der2<2>(0), support),
                1e-15);
}

TEST(Approximations, RBFBasis2DAtZero) {
    double s = 1.23;
    Gaussian<double> g(s);
    RBFBasis<Gaussian<double>, Vec2d> basis(2, g);
    Range<Vec2d> support = {{-0.77, 0.64}, {0.49, -2.1}};
    Vec2d p(0.0, 0.0);
    EXPECT_NEAR(basis.eval(0, p, support), basis.evalAt0(0, support), 1e-15);
    EXPECT_NEAR(basis.eval(1, p, support), basis.evalAt0(1, support), 1e-15);

    EXPECT_NEAR(basis.evalOp(0, p, Der1<2>(0), support),
                basis.evalOpAt0(0, Der1<2>(0), support, 1), 1e-14);
    EXPECT_NEAR(basis.evalOp(1, p, Der1<2>(1), support),
                basis.evalOpAt0(1, Der1<2>(1), support, 1), 1e-14);
    EXPECT_NEAR(basis.evalOp(0, p, Der2<2>(0, 1), support),
                basis.evalOpAt0(0, Der2<2>(0, 1), support, 1), 1e-14);
    EXPECT_NEAR(basis.evalOp(1, p, Der2<2>(0), support),
                basis.evalOpAt0(1, Der2<2>(0), support, 1), 1e-12);
    EXPECT_NEAR(basis.evalOp(0, p, Der2<2>(1), support),
                basis.evalOpAt0(0, Der2<2>(1), support, 1), 1e-12);
    EXPECT_NEAR(basis.evalOp(0, p, Lap<2>(), support),
                basis.evalOpAt0(0, Lap<2>(), support, 1), 1e-12);
}

TEST(Approximations, RBFBasisUsageExample) {
    /// [RBF basis usage example]
    Gaussian<double> g(1.3);
    RBFBasis<Gaussian<double>, Vec2d> basis(3, g);  // construct
    std::cout << basis << std::endl;
    Range<Vec2d> support = {{0.1, 2.3}, {0.2, 2.2}, {0.0, 2.25}};
    double val = basis.eval(1, 0.05, support);  // evaluate
    double der = basis.evalOp(1, 0.05, Lap<2>(), support);  // evaluate derivative
    /// [RBF basis usage example]
    EXPECT_DOUBLE_EQ(0.0640224990298612, val);
    EXPECT_DOUBLE_EQ(0.2649577880791847, der);
}

}  // namespace mm
