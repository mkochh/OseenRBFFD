#include <medusa/bits/operators/computeShapes.hpp>
#include <medusa/bits/operators/UniformShapeStorage.hpp>

#include <medusa/bits/domains/FindClosest.hpp>
#include <medusa/bits/operators/RaggedShapeStorage.hpp>
#include <medusa/bits/operators/ShapeStorage.hpp>
#include <medusa/bits/utils/stdtypesutils.hpp>
#include <medusa/bits/domains/DomainDiscretization.hpp>
#include <medusa/bits/domains/BoxShape.hpp>
#include <medusa/bits/approximations/WLS.hpp>
#include <medusa/bits/approximations/WeightFunction.hpp>
#include <medusa/bits/approximations/JacobiSVDWrapper_fwd.hpp>
#include <medusa/bits/approximations/Monomials.hpp>
#include <medusa/bits/approximations/ScaleFunction.hpp>
#include <medusa/bits/approximations/Operators_fwd.hpp>
#include <medusa/bits/utils/Timer.hpp>

#include "gtest/gtest.h"

namespace mm {


TEST(Operators, computeShapes) {
    /// [computeShapes usage example]
    BoxShape<Vec2d> box(0.0, 1.0);
    DomainDiscretization<Vec2d> d = box.discretizeWithStep(0.1);
    d.findSupport(FindClosest(9));
    WLS<Monomials<Vec2d>> approx(2);
    auto shapes = d.computeShapes<sh::lap|sh::d1>(approx);  // implicit call via domain hook

    // direct call which overwrites shapes for Laplacian for internal nodes in given storage
    sh::operator_tuple<sh::lap, 2>::type operators;
    computeShapes(d, approx, d.interior(), operators, &shapes);
    /// [computeShapes usage example]

    int node = 12;
    approx.compute(d.pos(node), d.supportNodes(node));
    Eigen::VectorXd sh = approx.getShape(Der1<2>(0));  // d/dx shape
    Eigen::VectorXd sh2 = shapes.d1(0, node);
    EXPECT_EQ(sh, sh2);

    node = 5;
    approx.compute(d.pos(node), d.supportNodes(node));
    sh = approx.getShape(Der1<2>(1));  // d/dy shape
    sh2 = shapes.d1(1, node);
    EXPECT_EQ(sh, sh2);

    node = 34;
    approx.compute(d.pos(node), d.supportNodes(node));
    sh = approx.getShape(Lap<2>());  // lap shape
    sh2 = shapes.laplace(node);
    EXPECT_EQ(sh, sh2);
}

TEST(Operators, computeShapesD2) {
    BoxShape<Vec2d> box(0.0, 1.0);
    DomainDiscretization<Vec2d> d = box.discretizeWithStep(0.1);
    d.findSupport(FindClosest(9));
    WLS<Monomials<Vec2d>> approx(2);
    auto shapes = d.computeShapes<sh::lap|sh::d2>(approx);

    int node = 22;
    approx.compute(d.pos(node), d.supportNodes(node));
    Eigen::VectorXd sh = approx.getShape(Der2<2>(0));  // d/dx^2 shape
    Eigen::VectorXd sh2 = shapes.d2(0, 0, node);
    EXPECT_EQ(sh, sh2);

    node = 0;
    approx.compute(d.pos(node), d.supportNodes(node));
    sh = approx.getShape(Der2<2>(1));  // d/dy^2 shape
    sh2 = shapes.d2(1, 1, node);
    EXPECT_EQ(sh, sh2);

    node = 46;
    approx.compute(d.pos(node), d.supportNodes(node));
    sh = approx.getShape(Der2<2>(0, 1));  // d/dy^2 shape
    sh2 = shapes.d2(0, 1, node);
    EXPECT_EQ(sh, sh2);

    node = 78;
    approx.compute(d.pos(node), d.supportNodes(node));
    sh = approx.getShape(Lap<2>());  // lap shape
    sh2 = shapes.laplace(node);
    EXPECT_EQ(sh, sh2);
}

/// [custom zero op]
// Dummy zero operator
struct Zero : Operator<Zero> {
    template <typename basis_t> typename basis_t::scalar_t applyAt0(
            const basis_t&, int, const std::vector<typename basis_t::vector_t>&,
            typename basis_t::scalar_t) const { return 0.0; }
};
/// [custom zero op]

TEST(Operators, ComputeShapesCustom) {
    BoxShape<Vec2d> box(0.0, 1.0);
    DomainDiscretization<Vec2d> d = box.discretizeWithStep(0.1);
    int n = 9;
    d.findSupport(FindClosest(n));
    WLS<Monomials<Vec2d>> approx(2);

    /// [custom op usage]
    std::tuple<Lap<2>, Der1s<2>, Zero> ops;
    auto shapes = d.computeShapes<decltype(ops)>(ops, approx);
    shapes.get<Zero>(3, 2);  // coefficient 4 of shape for zero in 3rd node.
    shapes.get<Der1s<2>>(0, 3, 4);  // coefficient 4 of shape for dx in 3rd node.
    /// [custom op usage]
    int N = d.size();
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < n; ++j) {
            EXPECT_EQ(0, shapes.get<Zero>(i, j));
        }
    }
}

TEST(Operators, computeShapesEmptySupport) {
    BoxShape<Vec2d> box(0.0, 1.0);
    DomainDiscretization<Vec2d> d = box.discretizeWithStep(0.1);
    WLS<Monomials<Vec2d>> approx(2);
    UniformShapeStorage<Vec2d, std::tuple<Lap<2>>> storage;
    sh::operator_tuple<sh::lap, 2>::type operators;
    EXPECT_DEATH(computeShapes(d, approx, d.interior(), operators, &storage),
            "Did you forget to find support before computing shapes?");
}

TEST(Operators, computeShapesCentralNotFirst) {
    double h = 0.1;
    Range<Vec2d> nodes = {{h, 0}, {-h, 0}, {0, h}, {0, -h}};
    DomainDiscretization<Vec2d> d = DomainDiscretization<Vec2d>(UnknownShape<Vec2d>{});
    int node = d.addInternalNode({0, 0}, 1);
    Range<int> support;
    for (auto p : nodes) {
        support.push_back(d.addInternalNode(p, 1));
    }
    d.support(node) = support;
    WLS<Monomials<Vec2d>> approx(1);
    auto shapes = d.computeShapes<sh::d1>(approx, {node});

    Eigen::VectorXd sh_x = shapes.d1(0, node);
    Eigen::VectorXd sh_y = shapes.d1(1, node);
    for (int i = 0; i < nodes.size(); ++i) {
        if (nodes[i][0] != 0) {
            EXPECT_NEAR(sh_x[i], signum(nodes[i][0]) * 0.5 / h, 1e-14);
            EXPECT_NEAR(sh_y[i], 0, 1e-14);
        } else {
            EXPECT_NEAR(sh_x[i], 0, 1e-14);
            EXPECT_NEAR(sh_y[i], signum(nodes[i][1]) * 0.5 / h, 1e-14);
        }
    }
}

}  // namespace mm
