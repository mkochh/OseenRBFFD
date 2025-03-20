#ifndef MEDUSA_BITS_IO_STL_FWD_HPP_
#define MEDUSA_BITS_IO_STL_FWD_HPP_

/**
 * @file
 * Support for reading stereolithography files.
 *
 * @example test/io/STL_test.cpp
 */

#include <cstdint>
#include <iosfwd>
#include <vector>
#include <cstring>

namespace mm {

/**
 * Support for reading stereolithography files: https://en.wikipedia.org/wiki/STL_(file_format)
 *
 * @note Reading headers from STL files is currently not supported.
 *
 * @snippet io/STL_test.cpp STL usage example
 * @ingroup io
 */
class STL {
  public:
    /// Holds one 3d Point in a STL file.
    struct Point { float x /** x coordinate */, y /** y coordinate */, z /** z coordinate */; };
    /// Holds one STL triangle, which consists of three points, a normal and an attribute.
    struct Triangle {
        Point normal /** triangle normal */, p1 /** first point */, p2 /** second point */,
              p3 /** third point */;
        uint16_t attribute;  ///< attribute of this triangle
    };

    /**
     * Read a binary STL file. Ascii files are currently not supported.
     * @param filename path to the file
     * @return A vector of triangles.
     */
    static std::vector<Triangle> read(const std::string& filename);
};

/// Print a STL point.
std::ostream& operator<<(std::ostream& os, const STL::Point& p);
/// Print a STL triangle.
std::ostream& operator<<(std::ostream& os, const STL::Triangle& v);


}  // namespace mm

#endif  // MEDUSA_BITS_IO_STL_FWD_HPP_
