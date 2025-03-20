#include <medusa/bits/io/STL.hpp>

/**
 * @file
 * Implementation of functions for working with stereolithography files.
 */

#include <fstream>
#include <stdexcept>
#include <medusa/bits/utils/assert.hpp>

namespace mm {

std::vector<STL::Triangle> STL::read(const std::string& filename) {
    std::ifstream file(filename, std::ios::in | std::ios::binary);
    assert_msg(file, "File '%s' could not be opened: %s", filename, strerror(errno));
    constexpr int HEADER_SIZE = 80;
    file.seekg(HEADER_SIZE);
    assert_msg(file, "Error reading header in file '%s': %s", filename, strerror(errno));
    uint32_t num_triangles;
    file.read(reinterpret_cast<char*>(&num_triangles), sizeof(uint32_t));
    assert_msg(file, "Error reading number of triangles in file '%s': %s", filename,
               strerror(errno));
    std::vector<Triangle> triangles(num_triangles);
    constexpr int TRIANGLE_BYTE_SIZE = 4 * 4 * 3;
    constexpr int ATTRIBUTE_BYTE_SIZE = 2;
    for (uint64_t i = 0; i < num_triangles; ++i) {
        file.read(reinterpret_cast<char*>(&triangles[i]), TRIANGLE_BYTE_SIZE);
        assert_msg(file, "Error reading triangle %d in file '%s': %s", i, filename,
                   strerror(errno));
        file.read(reinterpret_cast<char*>(&triangles[i].attribute), ATTRIBUTE_BYTE_SIZE);
        assert_msg(file, "Error reading attribute %d in file '%s': %s", i, filename,
                   strerror(errno));
    }
    return triangles;
}

std::ostream& operator<<(std::ostream& os, const STL::Point& p) {
    return os << "P(" << p.x << ", " << p.y << ", " << p.z << ")";
}
std::ostream& operator<<(std::ostream& os, const STL::Triangle& v) {
    return os << "T(" << v.normal << ", " << v.p1 << ", " << v.p2 << ", " << v.p3 << ")";
}

}  // namespace mm
