#include <medusa/bits/types/ScalarField.hpp>

#include "gtest/gtest.h"

namespace mm {

TEST(Types, ScalarFieldAssgn) {
    ScalarFieldd field(4);
    field = -3.0;
    EXPECT_EQ(4, field.size());
    EXPECT_EQ(-3, field[0]);
    EXPECT_EQ(-3, field[1]);
    EXPECT_EQ(-3, field[2]);
    EXPECT_EQ(-3, field[3]);

    indexes_t idx = {1, 2, 0, 0};
    field[idx] = 5.0;
    EXPECT_EQ(4, field.size());
    EXPECT_EQ(5.0, field[0]);
    EXPECT_EQ(5.0, field[1]);
    EXPECT_EQ(5.0, field[2]);
    EXPECT_EQ(-3.0, field[3]);

    ScalarFieldd field2 = field;
    field[{2, 1, 0}] = field2[{0, 2, 3}];
    EXPECT_EQ(4, field.size());
    EXPECT_EQ(-3.0, field[0]);
    EXPECT_EQ(5.0, field[1]);
    EXPECT_EQ(5.0, field[2]);
    EXPECT_EQ(-3.0, field[3]);
}

TEST(Types, ScalarFieldTraits) {
    static_assert(std::is_same<double, typename scalar_type<Eigen::Matrix2d>::type>::value, "");
    Eigen::MatrixX2f M;
    static_assert(std::is_same<float, typename scalar_type<
            decltype(M.row(0))>::type>::value, "");
    static_assert(std::is_same<float, typename scalar_type<
            decltype(M.block(0, 0, 1, 1))>::type>::value, "");
    static_assert(std::is_same<float, typename scalar_type<
            decltype(M({1, 2, 3}, {0, 1}))>::type>::value, "");
}

TEST(Types, DISABLED_ScalarFieldUsageExample) {
    /// [Scalar field usage example]
    ScalarFieldd field(4);
    // set to scalar and multi index access
    field = -3.0;
    field[{0, 2, 3}] = 5.0;
    // usual mathematical operations
    field *= 2.5;
    field += field;
    std::cout << field << std::endl;
    /// [Scalar field usage example]
}

}  // namespace mm
