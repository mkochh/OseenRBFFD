#include "medusa/bits/io/CSV.hpp"
#include "medusa/bits/io/CSV_Eigen.hpp"

#include "gtest/gtest.h"

#include <vector>
#include <Eigen/Core>

namespace mm {

TEST(IO, CSVRW1d) {
    std::vector<double> v = {1, 2, 3, -3.41203498123234};
    CSV::write("test/testdata/v.csv", v);

    std::vector<double> vread = CSV::read("test/testdata/v.csv");
    EXPECT_EQ(v, vread);
}

TEST(IO, CSVRW2d) {
    std::vector<std::vector<double>> M = {{1,    2},
                                          {3,    -3.4},
                                          {-4.5, 4.5}};
    CSV::write2d("test/testdata/M.csv", M, '!');
    auto Mread = CSV::read2d("test/testdata/M.csv", '!');
    EXPECT_EQ(M, Mread);
}

TEST(IO, CSVEigen1D) {
    Eigen::VectorXd v(10);
    v << 1, 3, 4, -2, 3, 5.6, 3.4, -2.3, 4, 5;
    CSV::write("test/testdata/v2.csv", v);
    auto vread = CSV::readEigen("test/testdata/v2.csv");

    EXPECT_TRUE((v - vread).isZero(0));
}

TEST(IO, CSVEigen2D) {
    Eigen::MatrixXd M(3, 2); M << 1, 2, 3, 4.3, 7.1, -2.3;
    CSV::writeEigen("test/testdata/M2.csv", M, ';');

    auto Mread = CSV::readEigen("test/testdata/M2.csv", ';');
    EXPECT_TRUE((M - Mread).isZero(0));
}

TEST(IO, CSVDeath) {
    EXPECT_DEATH(CSV::read2d("test/testdata/death1.csv"), "Failed parsing .* to double");
    EXPECT_DEATH(CSV::read2d("test/testdata/death2.csv"),
                 "Not all rows in CSV file have the same number of elements");
    EXPECT_DEATH(CSV::read2d("test/testdata/death3.csv"), "Error opening CSV file");
}

TEST(IO, CSVUsageExample) {
    /// [CSV usage example]
    std::vector<double> v = {1, 2.5, -4.5, 12.34};
    CSV::write("test/testdata/example.csv", v);

    std::vector<double> v2 = CSV::read("test/testdata/example.csv");
    // v is the same as v2

    // For the following, the CSV_Eigen.hpp header needs to be included.
    Eigen::MatrixXd M = CSV::readEigen("test/testdata/example.csv");
    CSV::writeEigen("test/testdata/example.csv", M.transpose()*M);
    /// [CSV usage example]
    EXPECT_EQ(v, v2);
}

}  // namespace mm
