#ifndef MEDUSA_BITS_DOMAINS_DOMAINDISCRETIZATION_FWD_HPP_
#define MEDUSA_BITS_DOMAINS_DOMAINDISCRETIZATION_FWD_HPP_

/**
 * @file
 * Declarations for domain discretizations.
 *
 * @example test/domains/DomainDiscretization_test.cpp
 */

#include <medusa/Config.hpp>
#include <medusa/bits/types/Range_fwd.hpp>
#include <medusa/bits/utils/memutils.hpp>
#include <medusa/bits/domains/DomainShape_fwd.hpp>
#include <medusa/bits/operators/shape_flags.hpp>
#include <iosfwd>
#include <string>

namespace mm {

template <typename vec_t, typename OpFamilies>
class RaggedShapeStorage;

/**
 * Class representing domain discretization along with an associated @ref DomainShape "shape".
 *
 * Domain discretization consists of points in the interior and on the boundary.
 * Each point has an associated type, with types of internal points being positive
 * and types of boundary points being negative. Boundary points also have an outer
 * unit normal normal associated with them.
 *
 * Domains can be @ref add "added" and @ref subtract "subtracted", nodes can be
 * @ref addInternalNode "added" and @ref removeNodes "removed" (removal currently
 * at O(N) cost). Support / stencil indices can be
 * @ref findSupport "stored" and @ref support "accessed".
 *
 * Additionally, domain interiors can be @ref GeneralFill "filled", @ref BasicRelax
 * "regularized" and used when @ref computeShapes "computing shapes".
 *
 * Usage example:
 * @snippet domains/DomainDiscretization_test.cpp Domain discretization usage example
 * @sa DomainShape, FindClosest
 * @ingroup domains
 */
template <class vec_t>
class DomainDiscretization {
  public:
    typedef vec_t vector_t;  ///< Vector data type used in computations.
    typedef typename vec_t::Scalar scalar_t;   ///< Scalar data type used in computation.
    /// Store dimension of the domain.
    enum { /** Dimensionality of the domain. */ dim = vec_t::dim };

  protected:
    Range<vec_t> positions_;  ///< Positions of internal discretization points.
    /**
     * Attribute used to store support points of each node. For each node `i`, `support[i]` is
     * list of indices of `i`'s support ordered by distance from `i`.
     * The first element in the `i`-th list is usually `i` followed by one or more others.
     **/
    Range<Range<int>> support_;
    /**
     * Storing types of nodes aligned with positions.
     * <b>Negative types are reserved for boundary nodes, positive for interior.
     * Zero type is not used.</b>
     * Example:
     * @code
     * positions_ = {p1, p2, p3, p4, p5};
     * types_ =     { 2, -1,  2,  1, -1};
     * @endcode
     * Nodes `p2` and `p5` are of type `-1`, `p1` and `p3` of type `2` and so on.
     */
    Range<int> types_;

    /**
     * Mapping index of a boundary node among all nodes, to its index among boundary nodes.
     * Indices of internal nodes are set to `-1`.
     * @code
     * positions_ =    {p1, p2, p3, p4, p5};
     * types_ =        { 2, -1,  2,  1, -1};
     * boundary_map_ = {-1,  0, -1, -1, 1};
     * @endcode
     * Nodes `p2` and `p5` are the first and the second boundary nodes.
     * There are no guarantees on the order of nodes.
     */
    Range<int> boundary_map_;

    /**
     * List of normals of boundary nodes.
     * @code
     * positions_ =    {p1, p2, p3, p4, p5};
     * types_ =        { 2, -1,  2,  1, -1};
     * boundary_map_ = {-1,  0, -1, -1, 1};
     * normals = {n1, n2};
     * @endcode
     * Normal `n1` is the normal of node `p2` and normal `n2` is the normal of `p5`.
     */
    Range<vec_t> normals_;

    /// Geometric shape of the domain.
    deep_copy_unique_ptr<DomainShape<vec_t>> shape_;

  public:
    /// Construct an empty discretization for a given shape.
    explicit DomainDiscretization(const DomainShape<vec_t>& shape);

    /**
     * Load discretization from a HD5 file.
     *
     * This function looks for datasets named `pos`, `types`, `normals` and `bmap`.
     * `pos` should be a `dim x N` dataset of doubles
     * `types` and `bmap` should be a `N` datasets of ints
     * `normals` should be a `dim x B` dataset of doubles, where `B` is the number of
     * boundary nodes. This format is produced by the HDF::writeDomain function.
     *
     * Usage example:
     * @snippet domains/DomainDiscretization_test.cpp load usage example
     *
     * @param file HDF object containing the discretization.
     * @param name Name of the folder where discretization is stored.
     *
     * @sa HDF::writeDomain
     */
    template <typename hdf5_file>
    static DomainDiscretization<vec_t> load(hdf5_file& file, const std::string& name);

    /// Returns positions of all nodes.
    const Range<vec_t>& positions() const { return positions_; }
    /// Returns the position of `i`-th node.
    const vec_t& pos(int i) const { return positions_[i]; }
    /// Returns writeable position of `i`-th node.
    vec_t& pos(int i) { return positions_[i]; }
    /// Returns `j`-th coordinate of the position of `i`-th node.
    scalar_t pos(int i, int j) const { return positions_[i][j]; }
    /// Returns support indices for all nodes.
    const Range<Range<int>>& supports() const { return support_; }
    /// Returns support indices for `i`-th node. Indices are ordered by distance to support nodes,
    /// with the first being the closest, usually the node itself.
    const Range<int>& support(int i) const { return support_[i]; }
    /// Returns writeable array of support indices for `i`-th node.
    Range<int>& support(int i) { return support_[i]; }
    /// Returns `j`-th support node of `i`-th node.
    int support(int i, int j) const { return support_[i][j]; }
    /// Returns positions of support nodes of `i`-th node.
    Range<vec_t> supportNodes(int i) const { return positions_[support_[i]]; }
    /// Returns position of `j`-th support node of `i`-th node.
    vec_t supportNode(int i, int j) const { return positions_[support_[i][j]]; }
    /// Returns Euclidean distance to the second support node.
    scalar_t dr(int i) const { return (positions_[i] - positions_[support_[i][1]]).norm(); }
    /// Returns size of `i`-th node support.
    int supportSize(int i) const { return support_[i].size(); }
    /// Returns a vector of support sizes for each node.
    Range<int> supportSizes() const;
    /// Returns types of all nodes.
    const Range<int>& types() const { return types_; }
    /// Returns mutable types of all nodes.
    Range<int>& types() { return types_; }
    /// Returns type of `i`-th node.
    int type(int i) const { return types_[i]; }
    /// Returns writeable type of `i`-th node.
    int& type(int i) { return types_[i]; }
    /// Returns geometric shape of the underlying domain.
    const DomainShape<vec_t>& shape() const { return *shape_; }

    /// Returns boundary map. @sa boundary_map_
    const Range<int>& bmap() const { return boundary_map_; }
    /**
     * Returns index of node `node` among only boundary nodes. The returned index is
     * in range `[0, boundary().size())` if `node` is a boundary node, and `-1` otherwise.
     */
    int bmap(int node) const { return boundary_map_[node]; }
    /// Returns normals of all boundary nodes.
    const Range<vec_t>& normals() const { return normals_; }
    /**
     * Returns outside unit normal of `i`-th node. The node must be a boundary node.
     * @throw Assertion fails if the noe is not a boundary node, i.e.\ `type(i) < 0` must hold.
     */
    const vec_t& normal(int i) const;
    /// Returns writable outside unit normal of `i`-th node. @sa normal
    vec_t& normal(int i);

    /// Returns indexes of all boundary nodes.
    indexes_t boundary() const { return types_ < 0; }
    /// Returns indexes of all internal nodes.
    indexes_t interior() const { return types_ > 0; }
    /// Returns indexes of all nodes, i.e.\ `{0, 1, ..., N-1}`.
    indexes_t all() const { return types_ != 0; }

    /// Returns `N`, the number of nodes in this discretization.
    int size() const { return positions_.size(); }

    /**
     * Adds a single interior node with specified type to this discretization.
     * @param point Coordinates of the node to add.
     * @param type Type of the node to add. Must be positive.
     * @return The index of the new node.
     * @sa addBoundaryNode
     */
    int addInternalNode(const vec_t& point, int type);

    /**
     * Adds a boundary node with given type and normal to the domain.
     * @param point Coordinates of the node to add.
     * @param type Type of the point, must be negative.
     * @param normal Outside unit normal to the boundary at point `point`.
     * @return The index of the new node.
     * @sa addInternalNode
     */
    int addBoundaryNode(const vec_t& point, int type, const vec_t& normal);

  private:
    /// Helper for adding points to the domain, which keeps the object consistent.
    int addNode(const vec_t& point, int type);
    /**
     * Remove a portion of this domain which is contained in `d` or is only `dx` away. The `dx`
     * is computed locally from the discretization of `d`. No nodes are added.
     *
     * Helper for @ref add and @ref subtract methods.
     *
     * @param d The discretization of the domain which specifies which part to chop off.
     * @param relevant Only the nodes with indices in `relevant` are used for computing the local
     * `dx`. This is mostly a performance improvement which allows to construct a smaller k-d tree.
     */
    void chopOff(const DomainDiscretization<vec_t>& d, const Range<int>& relevant);

  public:
    /**
     * Add nodes from another discretization to this discretization. The shape of the domain is not
     * changed. Returns indices of newly added nodes.
     */
    indexes_t addNodes(const DomainDiscretization<vec_t>& d);

    /// Changes node `i` to boundary point with given `type` and `normal`.
    void changeToBoundary(int i, const vec_t& point, int type, const vec_t& normal);

    /// Changes node `i` to interior point with given `type`.
    void changeToInterior(int i, const vec_t& point, int type);

    /// Overload of @ref addGhostNodes with constant `h` and adds ghost nodes to all boundary nodes.
    Range<int> addGhostNodes(scalar_t h, int type = 0);
    /// Overload of @ref addGhostNodes which adds ghost nodes to all boundary nodes.
    template <typename func_t>
    Range<int> addGhostNodes(func_t h, int type = 0);
    /// Overload of @ref addGhostNodes with constant `h`.
    Range<int> addGhostNodes(scalar_t h, int type, const indexes_t& indexes);
    /**
     * Add ghost or fictitious nodes to some boundary nodes.
     * @param h Nodal spacing function. The ghost node associated with some boundary node `p` and
     * with normal `n` is added at position `p + h(p)*n`.
     * @param type Type of the added ghost node. The usual choice is type `0`, making them neither
     * internal nor boundary. This means that shapes in these nodes will not be computed,
     * nor will their supports be found, however, they can appear as support nodes of other
     * domain nodes. Note that this means that their support size is 0, which has to be taken into
     * account when preallocating matrix sizes for implicit solving.
     * If a positive number is chosen for `type`, the nodes are added as internal nodes and
     * appear in the @ref interior index list. Similarly, if `type` is negative, they are added
     * as boundary nodes with normal `n` and appear in the @ref boundary index list.
     * @param indexes A list of indexes for which to add boundary nodes. If not given, all boundary
     * indices as given by @ref boundary are assumed.
     * @return A list of indices of the size of the original domain, mapping each node index
     * to the index of its corresponding ghost node, or -1, if it has no corresponding ghost node.
     * Precisely the slots with indices in `indexes` have an entry not equal to `-1`.
     *
     * This method can be useful when solving PDEs implicitly with Neumann boundary conditions,
     * see the [Ghost nodes](http://e6.ijs.si/medusa/wiki/index.php/Ghost_nodes) example
     * and the snippet below:
     * @snippet end2end/poisson_ghost_nodes.cpp Ghost nodes
     *
     * @warning Support nodes are not computed for ghost nodes, but ghost nodes are included
     * in other nodes' supports. This means that support size in ghost nodes will be zero,
     * also when returned by ShapeStorage::supportSizes, however, the matrix rows for ghost
     * nodes must still be preallocated.
     *
     * @sa addInternalNode, addBoundaryNode
     */
    template <typename func_t>
    Range<int> addGhostNodes(func_t h, int type, const indexes_t& indexes);

    /**
     * Shuffles node order according to given permutation. This also modifies the support indices
     * and all other necessary data.
     * For given list `I = {4, 3, 2, 0, 1}` and original positions `P = {a, b, c, d, e}`, the new
     * list of positions would be `{e, d, c, a, b}`, i.e. `P = P(I)`. This means that
     * `old_domain.pos(p[i]) == new_domain.pos(i)` holds for all nodes.
     * @param permutation A list of `n` integers containing each number from `0` to `n-1` exactly
     * once.
     * @throw Assertion fails if `permutation` is not a permutation.
     */
    void shuffleNodes(const indexes_t& permutation);

    /**
     * Reorders nodes according to the compare function. This can be useful to make the matrix
     * of the implicit system have a nicer structure. Sorts lexicographically (i.e. by x-coordinate
     * first) by default.
     * @param cmp A compare function for elements of type `std::pair<vec_t, int>` as passed
     * to `std::sort`. The first value in the pair is the coordinate of the point and the second
     * is its index.
     * @return The permutation used to shuffle the nodes.
     * @sa shuffleNodes
     */
    template <typename compare_t = std::less<std::pair<vec_t, int>>>
    Range<int> reorderNodes(const compare_t& cmp = compare_t());

    /**
     * Merges given discretization to the current discretization.
     * The shape of current domain changes to become the union of `shape()` and `d.shape()`.
     * The discretization `d` is added to `*this`. Types of nodes are preserved. This makes it
     * possible to set custom types before the `+=` operation and refer to those types in the
     * union as well.
     *
     * Example:
     * @snippet domains/DomainDiscretization_test.cpp add example
     *
     * More specifically, when performing `d1 += d2`, first, boundary nodes in `d2` that are
     * contained in `d1` are removed. Then, all nodes from `d1` that are included in
     * `d2` are removed. Additionally, nodes from `d1` are removed if  they are too close to
     * any of the nodes in `d2`. The threshold `dx` for point `p` in `d1` is computed by finding
     * the closest neighbour to `p` among the remaining nodes of `d1` and then finding its
     * closest neighbour among those same nodes. A factor of the distance between these last
     * two nodes is taken as `dx`. The remaining nodes in d2 are then added to d1.
     *
     * See the images below for visualization of the add operation.
     * @image html add_bnd.png
     * @image latex add_bnd.png
     * @image html add_full.png
     * @image latex add_full.png
     *
     * @sa ShapeUnion, subtract
     */
    DomainDiscretization& add(const DomainDiscretization& d);
    /// Operator version of DomainDiscretization::add. @sa add
    DomainDiscretization& operator+=(const DomainDiscretization& d) { return add(d); }

    /**
     * Subtracts given discretized domain from `*this`.
     *
     * The shape of current domain changes to become the difference of `shape()` and `d.shape()`.
     * The discretization nodes that are now outside of domain are removed
     * and the appropriate boundary nodes from `d` are added.  Types of nodes are preserved.
     * This makes it possible to set custom types before the `-=` operation and refer to those
     * types in the difference as well.
     *
     * Example:
     * @snippet domains/DomainDiscretization_test.cpp subtract example
     *
     * More specifically, when performing `d1 -= d2`, nodes in d1 that are contained
     * in `d2` or less than `dx` away from `d2` are removed. The local spacing `dx` for a point `p`
     * in `d1` is computed by first finding the nearest boundary node in `d2` and then taking a
     * fraction of the distance to its nearest boundary node.
     * Boundary nodes from `d2` that are contained in `d1` are also added, with inverted normals.
     *
     * See the images below for visualization of the subtract operation.
     * @image html sub_bnd.png
     * @image latex sub_bnd.png
     * @image html sub_full.png
     * @image latex sub_full.png
     *
     * @sa ShapeDifference, add
     */
    DomainDiscretization& subtract(const DomainDiscretization& d);
    /// Operator version of DomainDiscretization::subtract. @sa subtract
    DomainDiscretization& operator-=(const DomainDiscretization& d) { return subtract(d); }

    /// Translate the domain by a given vector `a`.
    DomainDiscretization& translate(const vec_t& a);

    /// Transform the domain by given orthogonal matrix `Q`.
    DomainDiscretization& rotate(const Eigen::Matrix<scalar_t, dim, dim>& Q);

    /// 2D version of @ref rotate accepting an angle.
    DomainDiscretization& rotate(scalar_t angle);

    /**
     * Checks if domain is in valid state. This includes that all points are contained in the
     * domain, that all boundary nodes have a normal, etc...
     */
    void assert_is_valid() const;

    /**
     * Clears all data about this discretization. The state of the object is as if it were
     * newly constructed using `shape()`.
     */
    void clear();

    /// Remove nodes with given indices. This removes their types and potential normals as well.
    void removeNodes(const Range<int>& to_remove);

    /// Removes all internal nodes.
    void removeInternalNodes() { removeNodes(types_ > 0); }

    /// Removes all boundary nodes.
    void removeBoundaryNodes() { removeNodes(types_ < 0); }

    /**
     * Returns `true` if `point` is inside the domain.
     * @throw Assertion fails if domain shape does not have a @ref contains method. In that
     * case, @ref discreteContains might be more appropriate.
     * @sa DomainShape::hasContains, DomainShape::contains, DomainDiscretization::discreteContains
     */
    bool contains(const vec_t& point) const;

    /**
     * A discrete version of @ref contains using given boundary discretization. A point `p` is
     * defined to be contained in the domain, iff the dot product `(p-q).n` is negative, where
     * `q` is the nearest boundary point and `n` is its outer normal.
     * This works sufficiently well, but can fail if boundary discretization is too sparse.
     * @param point The point being tested.
     * @param search Nearest neighbour search structure containing boundary points in the same order
     * as defined by `boundary_map_`. Such suitable structure is returned by
     * @ref makeDiscreteContainsStructure.
     *
     * Usage example:
     * @snippet DomainDiscretization_test.cpp discreteContains usage example
     * @sa contains, makeDiscreteContainsStructure
     */
    template <typename search_structure_t>
    bool discreteContains(const vec_t& point, const search_structure_t& search) const;

    /**
     * Returns `true` if `point` is inside the domain. If possible, this function uses
     * analytical DomainShape::contains, otherwise, it falls back to @ref discreteContains.
     * For parameter explanations see @ref discreteContains.
     * @sa contains, discreteContains, DomainShape::hasContains, DomainShape::contains
     */
    template <typename search_structure_t>
    bool contains(const vec_t& point, const search_structure_t& search) const;

    /**
     * Fills the search structure `search` with boundary points
     * for use with `contains()` or `discreteContains()`.
     * @param search A search structure to be used. Its contents are overwritten.
     * @sa discreteContains
     */
    template <typename search_structure_t>
    void makeDiscreteContainsStructure(search_structure_t& search) const;

    // TODO(jureslak): optimize order of nodes for matrix fill  -- sort by x coordinate or sth...
    // TODO(jureslak): symmetrize support?

    // ENGINE PLUGINS
    /// Define a method `Name` that calls its first argument.
    #define DOMAIN_PLUGIN(Name) \
        template<typename callable_t, typename... Args> \
        void Name(callable_t& callable, Args&&... args) { \
            return callable(*this, std::forward<Args>(args)...); \
        }

    /// Const version of DOMAIN_PLUGIN
    #define DOMAIN_PLUGIN_CONST(Name) \
        template<typename callable_t, typename... Args> \
        void Name(const callable_t& callable, Args&&... args) { \
            callable(*this, std::forward<Args>(args)...); \
        }

    /// Enables more readable calls to support engines.
    DOMAIN_PLUGIN(findSupport)
    /// Const version of DomainDiscretization::findSupport.
    DOMAIN_PLUGIN_CONST(findSupport)
    /// Enables more readable calls to fill engines.
    DOMAIN_PLUGIN(fill)
    /// Const version of DomainDiscretization::fill.
    DOMAIN_PLUGIN_CONST(fill)
    /// Enables more readable calls to relax engines.
    DOMAIN_PLUGIN(relax)
    /// Const version of DomainDiscretization::relax.
    DOMAIN_PLUGIN_CONST(relax)


    /// Compute shapes, specified with shape flags, for this domain with given approximation
    /// for given indexes. For complete documentation, se mm::computeShapes function.
    template <sh::shape_flags mask = sh::all, typename approx_t>
    RaggedShapeStorage<vec_t, typename sh::operator_tuple<mask, vec_t::dim>::type> computeShapes(
            approx_t approx, const indexes_t& indexes = {}) const;
    /// Compute shapes for default constructable families. @sa computeShapes
    template <typename OperatorTuple, typename approx_t>
    RaggedShapeStorage<vec_t, OperatorTuple> computeShapes(
            approx_t approx, const indexes_t& indexes = {}) const;
    /// Compute shapes for non-default constructable families. @sa computeShapes
    template <typename OperatorTuple, typename approx_t>
    RaggedShapeStorage<vec_t, OperatorTuple> computeShapes(OperatorTuple operators,
            approx_t approx, const indexes_t& indexes = {}) const;
    /// Output basic info about given domain. @sa computeShapes
    template <class V>
    friend std::ostream& operator<<(std::ostream& os, const DomainDiscretization<V>& d);
  private:
    /// Outputs a simple report about out domain, like number of nodes, support size.
    std::ostream& output_report(std::ostream& os = std::cout) const;
};

/// Output basic info about given domain.
template <class vec_t>
std::ostream& operator<<(std::ostream& os, const DomainDiscretization<vec_t>& d) {
    os << "Domain:\n";
    return d.output_report(os);
}

}  // namespace mm

#endif  // MEDUSA_BITS_DOMAINS_DOMAINDISCRETIZATION_FWD_HPP_
