#ifndef MEDUSA_BITS_DOMAINS_POLYHEDRONSHAPE_FWD_HPP_
#define MEDUSA_BITS_DOMAINS_POLYHEDRONSHAPE_FWD_HPP_

/**
 * @file
 * Implementation of shape representing polyhedrons.
 *
 * @example test/domains/PolyhedronShape_test.cpp
 */

#include <medusa/Config.hpp>
#include <medusa/bits/domains/discretization_helpers.hpp>
#include <medusa/bits/domains/discretization_helpers_advanced.hpp>
#include <medusa/bits/domains/GeneralFill.hpp>
#include <medusa/bits/domains/DomainShape.hpp>
#include <medusa/bits/domains/PolytopeShape.hpp>
#include <medusa/bits/types/Vec.hpp>

#include <CGAL/AABB_face_graph_triangle_primitive.h>

#ifdef USE_CGAL6
#include <CGAL/AABB_traits_3.h>
#endif

#ifdef USE_CGAL_LESS_THAN_6
#include <CGAL/AABB_traits.h>
#endif

#include <CGAL/AABB_tree.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

#include <fstream>

namespace mm {

/**
 * A polyhedron represented by a closed triangular mesh. This class uses CGAL
 * (https://www.cgal.org/) and its `Surface_mesh` class to do the heavy lifting.
 * @note CGAL (libcgal-dev) must be installed on your system and available in include and link
 * paths for this class to work.
 * @warning This header must be explicitly included to bet the 3D PolyhedronShape functionality. It is
 * not included by default due to its dependency on CGAL.
 *
 * @tparam vec_t This shape class is used only for 3D domains.
 *
 * Usage example:
 * @snippet domains/PolyhedronShape_test.cpp PolyhedronShape usage example
 * @ingroup domains
 */
template <typename vec_t>
class PolyhedronShape : public DomainShape<vec_t> {
    static_assert(vec_t::dim == 3, "Only available in 3 dimensions.");
    using base_t = DomainShape<vec_t>;  ///< Base class type.

    /* CGAL types */
    typedef CGAL::Simple_cartesian<double> K;  ///< CGAL geometry kernel.
    typedef K::Point_3 Point;  ///< CGAL point type.
    typedef K::Vector_3 Vector;  ///< CGAL vector type.
    typedef CGAL::Surface_mesh<Point> Mesh;  ///< CGAL surface mesh type.
    typedef Mesh::Vertex_index VI;  ///< CGAL index type for a point on a surface mesh.
    typedef Mesh::Face_index FI;  ///< CGAL index type for a face on a surface mesh.
    /// CGAL type for storing triangular faces of surface meshes in AABB tree.
    typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Primitive;
    /// CGAL type for axis-aligned bounding-box trees, used for testing point inclusions.
    #ifdef USE_CGAL6
    typedef CGAL::AABB_tree<CGAL::AABB_traits_3<K, Primitive>> Tree;
    #endif

    #ifdef USE_CGAL_LESS_THAN_6
    typedef CGAL::AABB_tree<CGAL::AABB_traits<K, Primitive>> Tree;
    #endif
    /// CGAL type of the function for testing point inclusion.
    typedef CGAL::Side_of_triangle_mesh<Mesh, K> PointInside;

    Mesh surface;  ///< Surface of the polyhedron represented as a triangular mesh.
    Tree tree;  ///< AABB tree of the faces.
    PointInside inside_tester;  ///< Function used for testing point inclusions.

    /// Initialize the tree and compute surface normals.
    void init() {
        tree.accelerate_distance_queries();
        auto vnormals = surface.add_property_map<VI, Vector>("v:normals", CGAL::NULL_VECTOR).first;
        auto fnormals = surface.add_property_map<FI, Vector>("f:normals", CGAL::NULL_VECTOR).first;
        CGAL::Polygon_mesh_processing::compute_normals(surface, vnormals, fnormals);
    }

  public:
    using typename base_t::scalar_t;
    using typename base_t::vector_t;
    using base_t::dim;
    using base_t::discretizeBoundaryWithDensity;
    using base_t::discretizeBoundaryWithStep;


    /// Construct a shape from a CGAL surface mesh. @sa fromOff.
    explicit PolyhedronShape(const Mesh& surface)
            : surface(surface),
              tree(surface.faces_begin(), surface.faces_end(), this->surface),
              inside_tester(tree) { init(); }

    /// Copy constructor (necessary to rebuild the tree).
    PolyhedronShape(const PolyhedronShape& other)
            : surface(other.surface),
              tree(surface.faces_begin(), surface.faces_end(), surface),
              inside_tester(tree) { init(); }

    /// Copy assignment (necessary to rebuild the tree).
    PolyhedronShape& operator=(const PolyhedronShape& other) {
        surface = other.surface;
        tree.clear();
        tree.insert(surface.faces_begin(), surface.faces_end(), surface);
        init();
        return *this;
    }

    /**
     * Read a triangular closed surface mesh describing a polyhedron from a `.off` file.
     * The face colors are read and stored as well and used to determine surface point types.
     * See the method @ref rgb2type for how RBG colors are converted to boundary types.
     * @param filename Path to the `.off` file.
     * @throws Assertion fails if file cannot be read, the surface is empty, not triangular or not
     * closed.
     */
    static PolyhedronShape fromOFF(const std::string& filename) {
        std::ifstream input(filename);
        assert_msg(input, "Failed opening file %s: %s", filename, strerror(errno));
        Mesh surface;
        input >> surface;
        assert_msg(!input.fail(), "Invalid OFF file - either parsing failed or the data "
                                  "does not represent a two-manifold.");
        assert_msg(input, "Failed reading file %s: %s", filename, strerror(errno));
        assert_msg(!surface.is_empty(), "Surface is empty.");
        assert_msg(CGAL::is_triangle_mesh(surface), "Surface must be a triangle mesh.");
        assert_msg(CGAL::is_closed(surface), "Surface must be closed (without holes).");
        return PolyhedronShape(surface);
    }

    bool contains(const vec_t& point) const override {
        return inside_tester(Point(point[0], point[1], point[2])) != CGAL::ON_UNBOUNDED_SIDE;
    }

    std::pair<vec_t, vec_t> bbox() const override {
        auto b = tree.bbox();
        return {{b.xmin(), b.ymin(), b.zmin()}, {b.xmax(), b.ymax(), b.zmax()}};
    }

    DomainDiscretization<vec_t> discretizeBoundaryWithDensity(
            const std::function<scalar_t(vec_t)>& dr, int type) const override;

    std::ostream& print(std::ostream& os) const override {
        return os << "PolyhedronShape with " << surface.num_vertices() << " vertices"
                  << " and " << surface.num_faces() << " faces.";
    }

    PolyhedronShape* clone() const override { return new PolyhedronShape(*this); }

    /// Convert a RBG color to an integer value that can be used as a boundary type.
    static int rgb2type(uint8_t r, uint8_t g, uint8_t b) {
        return -((static_cast<int>(r) << 16) + (static_cast<int>(g) << 8) + b) - 1;
    }

  private:
    /// Convert a CGAL vertex index to a Medusa vector.
    vec_t getPoint(VI idx) const {
        auto p = surface.point(idx);
        return {p.x(), p.y(), p.z()};
    }
};

// This must be hidden during a doxygen run due to a bug in processing const declarations with
// differently named template parameters.
// See https://github.com/doxygen/doxygen/issues/8178
#ifndef DOXYGEN_RUN
/**
 * Specialization for 3D, implementing polyhedrons.
 * @sa PolyhedronShape
 * @ingroup domains
 */
template <typename Scalar>
class PolytopeShape<Vec<Scalar, 3>> : public PolyhedronShape<Vec<Scalar, 3>> {
    using base_t = PolyhedronShape<Vec<Scalar, 3>>;
    using typename base_t::PolyhedronShape;
  public:
    PolytopeShape(const base_t& shape) : base_t(shape) {}
};
#endif


}  // namespace mm

#endif  // MEDUSA_BITS_DOMAINS_POLYHEDRONSHAPE_FWD_HPP_
