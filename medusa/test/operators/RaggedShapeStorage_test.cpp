#include <medusa/bits/types/Vec.hpp>
#include "medusa/bits/operators/RaggedShapeStorage.hpp"

#include "gtest/gtest.h"

namespace mm {

TEST(Operators, RaggedShapeStorageSuppSize) {
    RaggedShapeStorage<Vec2d, std::tuple<Lap<2>>> storage;
    Range<int> ss = {5, 2, 1};
    storage.resize(ss);

    Range<int> expected = ss;
    EXPECT_EQ(expected, storage.supportSizes());

    expected = {10, 4, 2, 10, 4, 2};
    EXPECT_EQ(expected, storage.supportSizesVec());

    RaggedShapeStorage<Vec3d, std::tuple<Lap<2>>> storage3;
    Range<int> ss3 = {5, 2, 1};
    storage3.resize(ss3);

    Range<int> expected3 = ss3;
    EXPECT_EQ(expected3, storage3.supportSizes());

    expected3 = {15, 6, 3, 15, 6, 3, 15, 6, 3};
    EXPECT_EQ(expected3, storage3.supportSizesVec());
}

TEST(Operators, RaggedShapeStorageSupp) {
    RaggedShapeStorage<Vec2d, std::tuple<Lap<2>>> storage;
    Range<int> ss = {1, 2, 5, 3, 1};
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

TEST(Operators, RaggedShapeStorageLap) {
    RaggedShapeStorage<Vec2d, std::tuple<Lap<2>>> storage;
    Range<int> ss = {1, 2, 5, 3, 1};
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

TEST(Operators, RaggedShapeStorageD1) {
    RaggedShapeStorage<Vec3d, std::tuple<Der1s<3>>> storage;
    Range<int> ss = {1, 2, 5, 3, 1};
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

TEST(Operators, RaggedShapeStorageD2) {
    RaggedShapeStorage<Vec3d, std::tuple<Der2s<3>>> storage;
    Range<int> ss = {1, 2, 5, 3, 1, 7, 12};
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


TEST(Operators, RaggedShapeStorageUsageExample) {
    /// [Ragged shape storage usage example]
    RaggedShapeStorage<Vec3d, std::tuple<Lap<3>, Der1s<3>>> storage;
    Range<int> sizes = {9, 13, 7};
    storage.resize(sizes);

    storage.size();  // 3
    Eigen::VectorXd lap(7);
    lap << 1.2, 3.4, 5.6, 7.8, 9.0, 1.2, 3.4;  // compute the shapes
    storage.setLaplace(2, lap);  // set lap as laplace shape for node 2.

    storage.laplace(2, 3);  // returns 7.8

    auto sh = storage.d1(1, 0);  // d/dy shape in node 0 (returns 0, because it is not set yet)
    // sh.size() is 9

    std::cout << storage << std::endl;
    /// [Ragged shape storage usage example]
    (void) (sh);  // otherwise unused
}

}  // namespace mm
