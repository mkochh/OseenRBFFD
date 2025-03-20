#include <medusa/bits/operators/UniformShapeStorage.hpp>
#include <medusa/bits/operators/shape_flags.hpp>
#include <medusa/bits/types/Vec.hpp>
#include "gtest/gtest.h"

namespace mm {

TEST(Operators, UniformShapeStorageSuppSize) {
    UniformShapeStorage<Vec2d, std::tuple<Lap<2>>> storage;
    Range<int> ss = {5, 5, 5};
    storage.resize(ss);

    Range<int> expected = ss;
    EXPECT_EQ(expected, storage.supportSizes());

    expected = {10, 10, 10, 10, 10, 10};
    EXPECT_EQ(expected, storage.supportSizesVec());
}

TEST(Operators, UniformShapeStorageSupp) {
    UniformShapeStorage<Vec2d, std::tuple<Lap<2>>> storage;
    Range<int> ss = {5, 5, 5};
    storage.resize(ss);

    for (int i = 0; i < ss.size(); ++i) {
        Range<int> supp(ss[i], i);
        storage.setSupport(i, supp);
        auto supp2 = storage.support(i);
        for (int j = 0; j < ss[i]; ++j) {
            EXPECT_EQ(supp[j], storage.support(i, j));
            EXPECT_EQ(supp[j], supp2[j]);
        }
    }
}

TEST(Operators, UniformShapeStorageLap) {
    UniformShapeStorage<Vec2d, std::tuple<Lap<2>>> storage;
    Range<int> ss = {5, 5, 5};
    storage.resize(ss);

    std::vector<Eigen::VectorXd> shapes;
    for (int i = 0; i < ss.size(); ++i) {
        Eigen::VectorXd sh(ss[i]);
        sh.setRandom();
        shapes.push_back(sh);
        storage.setLaplace(i, sh);
    }
    for (int i = 0; i < ss.size(); ++i) {
        Eigen::VectorXd sh2 = storage.laplace(i);
        EXPECT_EQ(shapes[i], sh2);

        for (int j = 0; j < ss[i]; ++j) {
            ASSERT_EQ(shapes[i][j], storage.laplace(i, j));
        }
    }
}

TEST(Operators, UniformShapeStorageD1) {
    UniformShapeStorage<Vec3d, std::tuple<Der1s<3>>> storage;
    Range<int> ss = {5, 5, 5, 5, 5, 5};
    storage.resize(ss);

    std::vector<std::vector<Eigen::VectorXd>> shapes(3);
    for (int d = 0; d < 3; ++d) {
        for (int i = 0; i < ss.size(); ++i) {
            Eigen::VectorXd sh(ss[i]);
            sh.setRandom();
            shapes[d].push_back(sh);
            storage.setD1(d, i, sh);
        }
    }
    for (int d = 0; d < 3; ++d) {
        for (int i = 0; i < ss.size(); ++i) {
            Eigen::VectorXd sh2 = storage.d1(d, i);
            const auto& sh = shapes[d][i];
            EXPECT_EQ(sh, sh2);
            for (int j = 0; j < ss[i]; ++j) {
                ASSERT_EQ(sh[j], storage.d1(d, i, j));
            }
        }
    }
}

TEST(Operators, UniformShapeStorageD2) {
    UniformShapeStorage<Vec3d, std::tuple<Der2s<3>>> storage;
    Range<int> ss = {6, 6, 6, 6, 6, 6, 6, 6};
    storage.resize(ss);

    std::vector<std::vector<std::vector<Eigen::VectorXd>>> shapes(3);
    for (int d = 0; d < 3; ++d) {
        for (int d2 = 0; d2 <= d; ++d2) {
            shapes[d].emplace_back();
            for (int i = 0; i < ss.size(); ++i) {
                Eigen::VectorXd sh(ss[i]);
                sh.setRandom();
                shapes[d].back().push_back(sh);
                storage.setD2(d2, d, i, sh);
            }
        }
    }

    for (int d = 0; d < 3; ++d) {
        for (int d2 = 0; d2 <= d; ++d2) {
            shapes[d].emplace_back();
            for (int i = 0; i < ss.size(); ++i) {
                const auto& sh = shapes[d][d2][i];
                Eigen::VectorXd sh2 = storage.d2(d2, d, i);
                EXPECT_EQ(sh, sh2);

                for (int j = 0; j < ss[i]; ++j) {
                    ASSERT_EQ(sh[j], storage.d2(d2, d, i, j));
                }
            }
        }
    }
}

TEST(Operators, DISABLED_UniformShapeStorageUsageExample) {
    /// [Uniform shape storage usage example]
    UniformShapeStorage<Vec3d, std::tuple<Lap<3>, Der1s<3>>> storage;
    Range<int> sizes = {5, 5, 5, 5, 5, 5};
    storage.resize(sizes);

    storage.size();  // 6
    Eigen::VectorXd lap(5);
    lap << 1.2, 3.4, 5.6, 7.8, 9.0;  // compute the shapes
    storage.setLaplace(2, lap);  // set lap as laplace shape for node 2.

    storage.laplace(2, 3);  // returns 7.8

    storage.d1(1, 3);  // d/dy shape in node 3 (returns 0, because it is not set yet)

    std::cout << storage << std::endl;
    /// [Uniform shape storage usage example]
}

}  // namespace mm
