#include "medusa/bits/domains/NURBSPatch.hpp"

#include "gtest/gtest.h"

namespace mm {

TEST(DomainEngines, NURBSPatchCurveEvaluate) {
    /// [NURBSPatch usage example]
    int p = 2;
    Range<double> w({1, sqrt(3) / 2, 1, sqrt(3) / 2, 1, sqrt(3) / 2, 1, sqrt(3) / 2, 1});
    Range<double> knots({0, 0, 0, PI / 2, PI / 2, PI, PI, 3 * PI / 2,
                         3 * PI / 2, 2 * PI, 2 * PI, 2 * PI});
    Range<Vec3d> cp({Vec3d(1, 0, 0), Vec3d(1, 1, 0), Vec3d(0, 1, 0), Vec3d(-1, 1, 0),
                     Vec3d(-1, 0, 0), Vec3d(-1, -1, 0), Vec3d(0, -1, 0), Vec3d(1, -1, 0),
                     Vec3d(1, 0, 0)});

    NURBSPatch<Vec3d, Vec1d> curve(cp, w, {knots}, {p});

    Vec3d pt = curve.evaluate({PI});
    /// [NURBSPatch usage example]
    EXPECT_LT((pt - Vec3d(-1, 0, 0)).norm(), 1e-15);

    pt = curve.evaluate({4.0});
    EXPECT_LT((pt - Vec3d(-0.680118, -0.779688, 0)).norm(), 1e-6);

    pt = curve.evaluate({0.25});
    EXPECT_LT((pt - Vec3d(0.973728, 0.266684, 0)).norm(), 1e-6);

    pt = curve.evaluate({0.0});
    EXPECT_LT((pt - Vec3d(1, 0, 0)).norm(), 1e-15);

    pt = curve.evaluate({2 * PI});
    EXPECT_LT((pt - Vec3d(1, 0, 0)).norm(), 1e-15);
}

TEST(DomainEngines, NURBSPatchCurveJacobian) {
    int p = 2;
    Range<double> w({1, sqrt(3) / 2, 1, sqrt(3) / 2, 1, sqrt(3) / 2, 1, sqrt(3) / 2, 1});
    Range<double> knots({0, 0, 0, PI / 2, PI / 2, PI, PI, 3 * PI / 2,
                         3 * PI / 2, 2 * PI, 2 * PI, 2 * PI});
    Range<Vec3d> cp({Vec3d(1, 0, 0), Vec3d(1, 1, 0), Vec3d(0, 1, 0), Vec3d(-1, 1, 0),
                     Vec3d(-1, 0, 0), Vec3d(-1, -1, 0), Vec3d(0, -1, 0), Vec3d(1, -1, 0),
                     Vec3d(1, 0, 0)});

    NURBSPatch<Vec3d, Vec1d> curve(cp, w, {knots}, {p});
    curve.computeDerivativeStructure();

    Vec3d pt = curve.jacobian({PI});
    EXPECT_LT((pt - Vec3d(0.0, -1.1026577908435842, 0.0)).norm(), 1e-15);

    pt = curve.jacobian({4.0});
    EXPECT_LT((pt - Vec3d(0.7398593320616775, -0.6222573803412518, 0.0)).norm(), 1e-6);

    pt = curve.jacobian({0.25});
    EXPECT_LT((pt - Vec3d(-0.2133476868599212, 1.0219703846566797, 0.0)).norm(), 1e-6);

    pt = curve.jacobian({0.0});
    EXPECT_LT((pt - Vec3d(0.0, 1.1026577908435842, 0.0)).norm(), 1e-15);

    pt = curve.jacobian({2 * PI});
    EXPECT_LT((pt - Vec3d(0.0, 1.1026577908435842, 0.0)).norm(), 1e-15);
}

TEST(DomainEngines, NURBSPatchCurveBoundary) {
    int p = 2;
    Range<double> w({1, sqrt(3) / 2, 1, sqrt(3) / 2, 1, sqrt(3) / 2, 1, sqrt(3) / 2, 1});
    Range<double> knots({0, 0, 0, PI / 2, PI / 2, PI, PI, 3 * PI / 2,
                         3 * PI / 2, 2 * PI, 2 * PI, 2 * PI});
    Range<Vec3d> cp({Vec3d(1, 0, 0), Vec3d(1, 1, 0), Vec3d(0, 1, 0), Vec3d(-1, 1, 0),
                     Vec3d(-1, 0, 0), Vec3d(-1, -1, 0), Vec3d(0, -1, 0), Vec3d(1, -1, 0),
                     Vec3d(1, 0, 0)});

    NURBSPatch<Vec3d, Vec1d> curve(cp, w, {knots}, {p});

    EXPECT_EQ(curve.getBoundaries().size(), 2);

    Vec1d pt = curve.getPatchParameterFromBoundaryParameter({}, 0, 0.25);
    EXPECT_LT((pt - Vec1d(PI / 2)).norm(), 1e-11);

    pt = curve.getPatchParameterFromBoundaryParameter({}, 1, 0.25);
    EXPECT_LT((pt - Vec1d(2 * PI - PI / 2)).norm(), 1e-11);
}

TEST(DomainEngines, NURBSPatchSurfaceEvaluate) {
    auto p = std::array<int, 2>{1, 2};
    auto knots = std::array<Range<double>, 2>{Range<double>({0, 0, 1, 1}),
                Range<double>({0, 0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1, 1})};
    Range<Range<Vec<double, 4>>> wcp(2);
    double x = 0.7071;

    wcp[0].push_back(Vec<double, 4>(1, 0, 0, 1));
    wcp[0].push_back(Vec<double, 4>(x, x, 0, x));
    wcp[0].push_back(Vec<double, 4>(0, 1, 0, 1));
    wcp[0].push_back(Vec<double, 4>(-x, x, 0, x));
    wcp[0].push_back(Vec<double, 4>(-1, 0, 0, 1));
    wcp[0].push_back(Vec<double, 4>(-x, -x, 0, x));
    wcp[0].push_back(Vec<double, 4>(0, -1, 0, 1));
    wcp[0].push_back(Vec<double, 4>(x, -x, 0, x));
    wcp[0].push_back(Vec<double, 4>(1, 0, 0, 1));

    wcp[1].push_back(Vec<double, 4>(1, 0, 1, 1));
    wcp[1].push_back(Vec<double, 4>(x, x, x, x));
    wcp[1].push_back(Vec<double, 4>(0, 1, 1, 1));
    wcp[1].push_back(Vec<double, 4>(-x, x, x, x));
    wcp[1].push_back(Vec<double, 4>(-1, 0, 1, 1));
    wcp[1].push_back(Vec<double, 4>(-x, -x, x, x));
    wcp[1].push_back(Vec<double, 4>(0, -1, 1, 1));
    wcp[1].push_back(Vec<double, 4>(x, -x, x, x));
    wcp[1].push_back(Vec<double, 4>(1, 0, 1, 1));

    NURBSPatch<Vec3d, Vec2d> surf(wcp, knots, p);

    Vec3d pt = surf.evaluate(Vec2d(0.5, 0.5));
    EXPECT_LT((pt - Vec3d(-1, 0, 0.5)).norm(), 1e-15);

    pt = surf.evaluate(Vec2d(0.543, 0.123));
    EXPECT_LT((pt - Vec3d(0.7164157130984911, 0.6976712959591224, 0.543)).norm(), 1e-15);

    pt = surf.evaluate(Vec2d(0.341, 0.796));
    EXPECT_LT((pt - Vec3d(0.26993127614132345, -0.9628790508533985, 0.341)).norm(), 1e-15);

    pt = surf.evaluate(Vec2d(0.0, 0.0));
    EXPECT_LT((pt - Vec3d(1.0, 0.0, 0.0)).norm(), 1e-15);

    pt = surf.evaluate(Vec2d(1.0, 1.0));
    EXPECT_LT((pt - Vec3d(1.0, 0.0, 1.0)).norm(), 1e-15);
}

TEST(DomainEngines, NURBSPatchSurfaceJacobian) {
    auto p = std::array<int, 2>{1, 2};
    auto knots = std::array<Range<double>, 2>{Range<double>({0, 0, 1, 1}),
            Range<double>({0, 0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1, 1})};
    Range<Range<Vec<double, 4>>> wcp(2);
    double x = 0.7071;

    wcp[0].push_back(Vec<double, 4>(1, 0, 0, 1));
    wcp[0].push_back(Vec<double, 4>(x, x, 0, x));
    wcp[0].push_back(Vec<double, 4>(0, 1, 0, 1));
    wcp[0].push_back(Vec<double, 4>(-x, x, 0, x));
    wcp[0].push_back(Vec<double, 4>(-1, 0, 0, 1));
    wcp[0].push_back(Vec<double, 4>(-x, -x, 0, x));
    wcp[0].push_back(Vec<double, 4>(0, -1, 0, 1));
    wcp[0].push_back(Vec<double, 4>(x, -x, 0, x));
    wcp[0].push_back(Vec<double, 4>(1, 0, 0, 1));

    wcp[1].push_back(Vec<double, 4>(1, 0, 1, 1));
    wcp[1].push_back(Vec<double, 4>(x, x, x, x));
    wcp[1].push_back(Vec<double, 4>(0, 1, 1, 1));
    wcp[1].push_back(Vec<double, 4>(-x, x, x, x));
    wcp[1].push_back(Vec<double, 4>(-1, 0, 1, 1));
    wcp[1].push_back(Vec<double, 4>(-x, -x, x, x));
    wcp[1].push_back(Vec<double, 4>(0, -1, 1, 1));
    wcp[1].push_back(Vec<double, 4>(x, -x, x, x));
    wcp[1].push_back(Vec<double, 4>(1, 0, 1, 1));

    NURBSPatch<Vec3d, Vec2d> surf(wcp, knots, p);
    surf.computeDerivativeStructure();

    Eigen::Matrix<double, 3, 2> jm = surf.jacobian(Vec2d(0.5, 0.5));
    Eigen::Matrix<double, 3, 2> actual_jm;
    actual_jm.col(0) << 0, 0, 1;
    actual_jm.col(1) << 0.0, -5.6568, 0.0;
    EXPECT_LT((jm - actual_jm).norm(), 1e-15);

    jm = surf.jacobian(Vec2d(0.543, 0.123));
    actual_jm.col(0) << 0, 0, 1;
    actual_jm.col(1) << -4.623582182148408, 4.747803096778865, 0.0;
    EXPECT_LT((jm - actual_jm).norm(), 1e-10);

    jm = surf.jacobian(Vec2d(0.341, 0.796));
    actual_jm.col(0) << 0, 0, 1;
    actual_jm.col(1) << 5.972112805677801, 1.6742281610437675, 0.0;
    EXPECT_LT((jm - actual_jm).norm(), 1e-10);

    jm = surf.jacobian(Vec2d(0.0, 0.0));
    actual_jm.col(0) << 0, 0, 1;
    actual_jm.col(1) << 0.0, 5.6568, 0.0;
    EXPECT_LT((jm - actual_jm).norm(), 1e-15);

    jm = surf.jacobian(Vec2d(1.0, 1.0));
    actual_jm.col(0) << 0, 0, 1;
    actual_jm.col(1) << 0.0, 5.6568, 0.0;
    EXPECT_LT((jm - actual_jm).norm(), 1e-15);
}

TEST(DomainEngines, NURBSPatchSurfaceBoundary) {
    auto p = std::array<int, 2>{1, 2};
    auto knots = std::array<Range<double>, 2>{Range<double>({0, 0, 1, 1}),
    Range<double>({0, 0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1, 1})};
    Range<Range<Vec<double, 4>>> wcp(2);
    double x = 0.7071;

    wcp[0].push_back(Vec<double, 4>(1, 0, 0, 1));
    wcp[0].push_back(Vec<double, 4>(x, x, 0, x));
    wcp[0].push_back(Vec<double, 4>(0, 1, 0, 1));
    wcp[0].push_back(Vec<double, 4>(-x, x, 0, x));
    wcp[0].push_back(Vec<double, 4>(-1, 0, 0, 1));
    wcp[0].push_back(Vec<double, 4>(-x, -x, 0, x));
    wcp[0].push_back(Vec<double, 4>(0, -1, 0, 1));
    wcp[0].push_back(Vec<double, 4>(x, -x, 0, x));
    wcp[0].push_back(Vec<double, 4>(1, 0, 0, 1));

    wcp[1].push_back(Vec<double, 4>(1, 0, 1, 1));
    wcp[1].push_back(Vec<double, 4>(x, x, x, x));
    wcp[1].push_back(Vec<double, 4>(0, 1, 1, 1));
    wcp[1].push_back(Vec<double, 4>(-x, x, x, x));
    wcp[1].push_back(Vec<double, 4>(-1, 0, 1, 1));
    wcp[1].push_back(Vec<double, 4>(-x, -x, x, x));
    wcp[1].push_back(Vec<double, 4>(0, -1, 1, 1));
    wcp[1].push_back(Vec<double, 4>(x, -x, x, x));
    wcp[1].push_back(Vec<double, 4>(1, 0, 1, 1));

    NURBSPatch<Vec3d, Vec2d> surf(wcp, knots, p);
    surf.computeDerivativeStructure();
    auto edges = surf.getBoundaries();

    double t = 0.3;
    for (int i = 0; i < 4; i++) {
        Vec3d pt_surf = surf.evaluate(surf.getPatchParameterFromBoundaryParameter(t, i));
        Vec3d pt_edge = edges[i].evaluate(t);
        EXPECT_LT((pt_surf - pt_edge).norm(), 1e-9);
    }

    t = 0.75;
    for (int i = 0; i < 4; i++) {
        Vec3d pt_surf = surf.evaluate(surf.getPatchParameterFromBoundaryParameter(t, i));
        Vec3d pt_edge = edges[i].evaluate(t);
        EXPECT_LT((pt_surf - pt_edge).norm(), 1e-9);
    }

    edges[0].computeDerivativeStructure();
    Vec3d der = edges[0].jacobian(t);
    auto jac = surf.jacobian(surf.getPatchParameterFromBoundaryParameter(t, 0));
    EXPECT_LT((der - jac.col(1)).norm(), 1e-9);

    edges[3].computeDerivativeStructure();
    der = edges[3].jacobian(t);
    jac = surf.jacobian(surf.getPatchParameterFromBoundaryParameter(t, 3));
    EXPECT_LT((der - jac.col(0)).norm(), 1e-9);
}

}  // namespace mm
