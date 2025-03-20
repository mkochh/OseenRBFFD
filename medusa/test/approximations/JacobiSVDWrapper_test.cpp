#include "medusa/bits/approximations/JacobiSVDWrapper.hpp"
#include <Eigen/SVD>

#include "gtest/gtest.h"

namespace mm {

TEST(Approximations, JacobiSVDWrapper) {
    /// [JacobiSVDWrapper usage example]
    Eigen::MatrixXd M(7, 9);
    M.setRandom();

    JacobiSVDWrapper<double> svd;
    svd.compute(M);

    // This code is the save as above, with explicit parameter specification.
    Eigen::JacobiSVD<Eigen::MatrixXd> svd2;
    svd2.compute(M, Eigen::ComputeThinU | Eigen::ComputeThinV);

    // Use svd.matrixU(), svd.solve(b) as usual.
    /// [JacobiSVDWrapper usage example]

    EXPECT_TRUE((svd.matrixU() - svd2.matrixU()).isZero(0));
    EXPECT_TRUE((svd.matrixV() - svd2.matrixV()).isZero(0));
}

}  // namespace mm
