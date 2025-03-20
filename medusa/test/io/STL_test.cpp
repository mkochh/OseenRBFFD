#include <medusa/bits/io/STL.hpp>
#include <medusa/bits/utils/print.hpp>
#include <medusa/bits/utils/assert.hpp>

#include "gtest/gtest.h"

namespace mm {

bool operator==(const STL::Point& a, const STL::Point& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z;
}

TEST(IO, STLReadSample) {
    /// [STL usage example]
    std::vector<STL::Triangle> stl = STL::read("test/testdata/tetrahedron.stl");
    std::stringstream ss; ss << stl;
    /// [STL usage example]

    std::string output = "["
         "T(P(0, 0, 1), P(0.5, 0.5, 0.5), P(1.5, 0.5, 0.5), P(0.5, 1.5, 0.5)), "
         "T(P(0, -1, 0), P(0.5, 0.5, 0.5), P(1.5, 0.5, 0.5), P(0.5, 0.5, 1.5)), "
         "T(P(1, 0, 0), P(0.5, 0.5, 0.5), P(0.5, 1.5, 0.5), P(0.5, 0.5, 1.5)), "
         "T(P(0.57735, 0.57735, 0.57735), P(1.5, 0.5, 0.5), P(0.5, 1.5, 0.5), P(0.5, 0.5, 1.5))]";
    EXPECT_EQ(output, ss.str());

    std::vector<STL::Point> vertices = {
            {0.5, 0.5, 0.5}, {1.5, 0.5, 0.5}, {0.5, 1.5, 0.5}, {0.5, 0.5, 1.5}};
    std::vector<std::vector<int>> faces = {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}};
    float f = 0.57735026918962576451f;
    std::vector<STL::Point> normals = {
            {0.0, 0.0, 1.0}, {0.0, -1.0, 0.0}, {1.0, 0.0, 0.0}, {f, f, f}};
    ASSERT_EQ(faces.size(), stl.size());
    for (int i = 0; i < static_cast<int>(faces.size()); ++i) {
        EXPECT_EQ(vertices[faces[i][0]], stl[i].p1);
        EXPECT_EQ(vertices[faces[i][1]], stl[i].p2);
        EXPECT_EQ(vertices[faces[i][2]], stl[i].p3);
        EXPECT_EQ(normals[i], stl[i].normal);
        EXPECT_EQ(i+1, stl[i].attribute);
    }
}

TEST(IO, STLreadUtah) {
    auto res = STL::read("test/testdata/Utah_teapot.stl");
    EXPECT_EQ(9438, res.size());
}

}  // namespace mm
