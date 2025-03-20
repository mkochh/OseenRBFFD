#ifndef MEDUSA_BITS_DOMAINS_STLSHAPE_FWD_HPP_
#define MEDUSA_BITS_DOMAINS_STLSHAPE_FWD_HPP_

#include "DomainShape_fwd.hpp"
#include <medusa/bits/io/STL_fwd.hpp>
#include <utility>

/**
 * @file
 * Declaration of STLShape class.
 *
 * @example test/domains/STLShape_test.cpp
 */

namespace mm {

/**
 * Class representing an object given by the STL file.
 *
 * This class offers basic support for STL files, with rudimentary surface discretization and
 * inexact geometry for the @ref contains method. For a more sophisticated class see the 3D version
 * of @ref PolygonShape, which uses the CGAL library that offers support for simple loading of
 * models from more file types, exact geometry for @ref contains and support for designating types
 * to faces.
 *
 * @snippet domains/STLShape_test.cpp STL shape usage example
 * @ingroup domains
 */
template <typename vec_t>
class STLShape : public DomainShape<vec_t> {
  public:
    using typename DomainShape<vec_t>::scalar_t;
    using DomainShape<vec_t>::dim;

  private:
    static_assert(dim == 3, "STL files are for 3D models only.");
    std::vector<vec_t> vertices_;  ///< 3d vertices.
    std::vector<std::array<int, 3>> faces_;  ///< Faces, described with three indices of vertices.
    std::vector<vec_t> normals_;  ///< Normals for each face.
    std::pair<vec_t, vec_t> bbox_;  ///< Bounding box.

  public:
    /// Create STLShape from triangles read from STL file.
    STLShape(const std::vector<STL::Triangle>& stl);

    /**
     * Contains method is not supported and is approximated with a discrete version, which can
     * be wrong.
     * @sa DomainDiscretization::discreteContains.
     */
    bool contains(const vec_t& /* point */) const override { return true; }

    bool hasContains() const override { return false; }

    std::pair<vec_t, vec_t> bbox() const override { return bbox_; }

    using DomainShape<vec_t>::discretizeBoundaryWithDensity;
    DomainDiscretization<vec_t> discretizeBoundaryWithDensity(
            const std::function<scalar_t(vec_t)>&, int) const override;

    std::ostream& print(std::ostream& ostream) const override {
        return ostream << "STL shape with " << vertices_.size() << " vertices and "
                       << faces_.size() << " faces.";
    }

    STLShape* clone() const override {
        return new STLShape<vec_t>(*this);
    }

  private:
    /// Convert STL point to Vec3d.
    static vec_t v(const STL::Point& p) {
        return vec_t(p.x, p.y, p.z);
    }
};

}  // namespace mm

#endif  // MEDUSA_BITS_DOMAINS_STLSHAPE_FWD_HPP_
