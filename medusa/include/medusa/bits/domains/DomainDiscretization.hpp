#ifndef MEDUSA_BITS_DOMAINS_DOMAINDISCRETIZATION_HPP_
#define MEDUSA_BITS_DOMAINS_DOMAINDISCRETIZATION_HPP_

/**
 * @file
 * Implementation of discretized domains.
 */

#include <medusa/bits/utils/assert.hpp>
#include <medusa/bits/spatial_search/KDTree.hpp>
#include <cmath>
#include "DomainDiscretization_fwd.hpp"
#include "DomainShape.hpp"
#include "ShapeUnion.hpp"
#include "ShapeDifference.hpp"
#include "UnknownShape.hpp"

namespace mm {

template <class vec_t>
Range<int> DomainDiscretization<vec_t>::supportSizes() const {
    int N = size();
    Range<int> sizes(N);
    for (int i = 0; i < N; ++i) {
        sizes[i] = supportSize(i);
    }
    return sizes;
}

template <class vec_t>
vec_t& DomainDiscretization<vec_t>::normal(int i) {
    assert_msg(0 <= i && i <= size(), "Index %d out of range [0, %d)", i, size());
    assert_msg(types_[i] < 0, "Node %d must be a boundary node, got type %d.", i, types_[i]);
    assert_msg(boundary_map_[i] != -1, "Node %d does not have a normal. Maybe you manually set"
            " supports and positions instead of using addInternalNode* methods?", i);
    return normals_[boundary_map_[i]];
}

template <class vec_t>
const vec_t& DomainDiscretization<vec_t>::normal(int i) const {
    assert_msg(0 <= i && i <= size(), "Index %d out of range [0, %d)", i, size());
    assert_msg(types_[i] < 0, "Node %d must be a boundary node, got type %d.", i, types_[i]);
    assert_msg(boundary_map_[i] != -1, "Node %d does not have a normal. Maybe you manually set"
            " supports and positions instead of using addInternalNode* methods?", i);
    return normals_[boundary_map_[i]];
}

template <class vec_t>
int DomainDiscretization<vec_t>::addInternalNode(const vec_t& point, int type) {
    assert_msg(type > 0, "This function is for adding internal points, but got type %d, which is "
                         "not positive. Use addBoundaryNode to add boundary nodes.", type);
    return addNode(point, type);
}

template <class vec_t>
int DomainDiscretization<vec_t>::addBoundaryNode(
        const vec_t& point, int type, const vec_t& normal) {
    assert_msg(type < 0, "Type of boundary points must be negative, got %d.", type);
    int idx = addNode(point, type);
    boundary_map_[idx] = normals_.size();
    normals_.push_back(normal);
    return idx;
}

template <class vec_t>
int DomainDiscretization<vec_t>::addNode(const vec_t& point, int type) {
    positions_.push_back(point);
    types_.push_back(type);
    support_.emplace_back();
    boundary_map_.push_back(-1);
    return positions_.size()-1;
}

template <class vec_t>
void DomainDiscretization<vec_t>::changeToBoundary(int i, const vec_t& point, int type,
                                                   const vec_t& normal)  {
    assert_msg(0 <= i && i < size(), "Index %d out of range [0, %d).", i, size());
    assert_msg(types_[i] >= 0, "Point %d is already a boundary point with type %d.", i,
               types_[i]);
    assert_msg(type < 0, "New type must be negative, got %d.", type);
    positions_[i] = point;
    types_[i] = type;
    boundary_map_[i] = normals_.size();
    normals_.push_back(normal);
}

template <class vec_t>
void DomainDiscretization<vec_t>::changeToInterior(int i, const vec_t& point, int type)  {
    assert_msg(0 <= i && i < size(), "Index %d out of range [0, %d).", i, size());
    assert_msg(types_[i] <= 0, "Point %d is already an interior point with type %d.", i, types_[i]);
    assert_msg(type > 0, "New type must be positive, got %d.", type);
    positions_[i] = point;
    types_[i] = type;
    normals_.remove({boundary_map_[i]});
    for (int idx = 0; idx < boundary_map_.size(); ++idx) {
        if (boundary_map_[idx] != -1 && boundary_map_[idx] > boundary_map_[i]) {
            boundary_map_[idx] -= 1;
        }
    }
    boundary_map_[i] = -1;
}

template <class vec_t>
indexes_t DomainDiscretization<vec_t>::addNodes(const DomainDiscretization<vec_t>& d) {
    indexes_t indexes(d.size());
    for (int i = 0; i < d.size(); ++i) {
        if (d.type(i) < 0) {
            indexes[i] = addBoundaryNode(d.pos(i), d.type(i), d.normal(i));
        } else {
            indexes[i] = addInternalNode(d.pos(i), d.type(i));
        }
    }
    return indexes;
}

template <class vec_t>
void DomainDiscretization<vec_t>::removeNodes(const Range<int>& to_remove)  {
    int n = size();
    Range<bool> remove_mask(n, false);
    remove_mask[to_remove] = true;
    Range<int> remove_normals;
    for (int i = 0; i < n; ++i) {
        if (boundary_map_[i] >= 0 && remove_mask[i]) {  // boundary node that has to be removed
            remove_normals.push_back(boundary_map_[i]);
        } else if (boundary_map_[i] >= 0) {
            boundary_map_[i] -= remove_normals.size();
        }
    }
    normals_.remove(remove_normals);
    boundary_map_.remove(to_remove);
    types_.remove(to_remove);
    positions_.remove(to_remove);
    if (support_.size() == n) support_.remove(to_remove);
}

template <class vec_t>
Range<int> DomainDiscretization<vec_t>::addGhostNodes(scalar_t h, int type) {
    return addGhostNodes([=](const vec_t&) { return h; }, type);
}

template <class vec_t>
Range<int> DomainDiscretization<vec_t>::addGhostNodes(
        scalar_t h, int type, const indexes_t& indexes) {
    return addGhostNodes([=](const vec_t&) { return h; }, type, indexes);
}

template <class vec_t>
template <typename func_t>
Range<int> DomainDiscretization<vec_t>::addGhostNodes(func_t h, int type) {
    return addGhostNodes(h, type, boundary());
}

template <class vec_t>
template <typename func_t>
Range<int> DomainDiscretization<vec_t>::addGhostNodes(func_t h, int type, const indexes_t& idx) {
    Range<int> gh(size(), -1);
    for (int i : idx) {
        assert_msg(0 <= i && i < size(), "Index %d out of range [0, %d)", i, size());
        assert_msg(types_[i] < 0, "Only boundary nodes can have associated ghost nodes, "
                                  "node %d has type %d.", i, types_[i]);
        scalar_t ch = h(pos(i));
        if (type > 0) gh[i] = addInternalNode(pos(i) + ch*normal(i), type);
        else if (type < 0) gh[i] = addBoundaryNode(pos(i) + ch*normal(i), type, normal(i));
        else if (type == 0) gh[i] = addNode(pos(i) + ch*normal(i), type);
    }
    return gh;
}

template <class vec_t>
void DomainDiscretization<vec_t>::clear() {
    positions_.clear();
    types_.clear();
    support_.clear();
    normals_.clear();
    boundary_map_.clear();
}

template <class vec_t>
void DomainDiscretization<vec_t>::assert_is_valid() const {
    int N = size();
    assert_msg(N == types_.size(), "Number of node type entries (%d) not equal to domain size "
            "(%d).", types_.size(), N);
    assert_msg(N == boundary_map_.size(), "Number of boundary map entries (%d) not equal to domain"
            " size (%d).", types_.size(), N);
    assert_msg(N == support_.size(), "Number of node supports (%d) not equal to domain size "
            "(%d).", types_.size(), N);
    assert_msg(N == positions_.size(), "Number of node positions (%d) not equal to domain size "
            "(%d).", types_.size(), N);
    int num_bnd = 0;
    for (int i = 0; i < N; ++i) {
        if (types_[i] < 0) {
            assert_msg(0 <= boundary_map_[i] && boundary_map_[i] < normals_.size(),
                       "Boundary map index %d of node %d must lie in interval [0, %d).",
                       boundary_map_[i], i, normals_.size());
            ++num_bnd;
        } else if (types_[i] > 0) {
            assert_msg(boundary_map_[i] == -1, "Expected boundary map of internal node %d "
                    "to be -1, got %d.", i, boundary_map_[i]);
        }
        assert_msg(shape_->contains(pos(i)), "Position %s at index %d is not contained "
                "in the domain.", pos(i), i);
    }
    assert_msg(num_bnd == normals_.size(), "Number of boundary nodes not the same as number "
            "of normals, got %d nodes with negative type, but %d normals.",
               num_bnd, normals_.size());
}

template <class vec_t>
std::ostream& DomainDiscretization<vec_t>::output_report(std::ostream& os) const  {
    int internal_count = (types_ > 0).size();
    int boundary_count = (types_ < 0).size();
    std::size_t mem = sizeof(*this);
    std::size_t mem_pos = mem_used(positions_);
    std::size_t mem_types = mem_used(types_);
    std::size_t mem_supp = mem_used(support_);
    int max_support_size = -1;
    for (int i = 0; i < support_.size(); ++i) {
        max_support_size = std::max(max_support_size, support_[i].size());
        mem_supp += mem_used(support_[i]);
    }
    std::size_t mem_bnd_map = mem_used(boundary_map_);
    std::size_t mem_normals = mem_used(normals_);
    std::size_t mem_total = mem + mem_pos + mem_types + mem_supp + mem_bnd_map + mem_normals;
    os << "    dimension: " << dim << '\n'
       << "    shape: " << *shape_ << '\n'
       << "    number of points: " << size() << " of which " << internal_count
       << " internal and " << boundary_count << " boundary\n"
       << "    max support size: ";
    if (max_support_size == -1) os << "support not found yet";
    else os << max_support_size;
    os << "\n    mem used total: " << mem2str(mem_total) << '\n'
       << "        pos & types:   " << mem2str(mem_pos + mem_types) << '\n'
       << "        supp:   " << mem2str(mem_supp) << '\n'
       << "        bnd & normals: " << mem2str(mem_bnd_map + mem_normals) << '\n';
    return os;
}

template <class vec_t>
DomainDiscretization<vec_t>::DomainDiscretization(const DomainShape<vec_t>& shape) :
        shape_(shape) {}

template <class vec_t>
void DomainDiscretization<vec_t>::chopOff(const DomainDiscretization& d,
                                          const Range<int>& relevant) {
    /* Chops off `d` from the current domain. Nodes in *this that are contained
     * in `d` or are less than `0.75*dx` away from `d` are removed. No nodes are added. */

    const bool tree_ok = relevant.size() >= 2;
    KDTree<vec_t> tree(d.positions_[relevant]);

    Range<int> to_remove;
    for (int i = 0; i < size(); ++i) {
        if (d.shape().contains(pos(i))) {
            to_remove.push_back(i);
            continue;
        }
        bool not_ok = false;
        if (tree_ok) {
            auto r = tree.query(pos(i));
            int idx = r.first[0];
            scalar_t dist = r.second[0];
            scalar_t dx = tree.query(d.pos(relevant[idx]), 2).second[1];
            not_ok = dist < 0.75*dx;
        }
        if (not_ok) {
            to_remove.push_back(i);
        }
    }
    removeNodes(std::move(to_remove));
}

template <class vec_t>
DomainDiscretization<vec_t>&
DomainDiscretization<vec_t>::add(const DomainDiscretization<vec_t>& d) {
    /*
     * Domain 1 is chopped off, and the relevant part of `d` is added. The added part is almost
     * the whole domain, without any boundary nodes that are contained in domain 1.
     */
    Range<int> to_add;
    for (int i = 0; i < d.size(); ++i) {
        if (d.type(i) > 0 || !shape().contains(d.pos(i))) {
            to_add.push_back(i);
        }
    }
    chopOff(d, to_add);
    for (int i : to_add) {
        if (d.type(i) < 0) {
            addBoundaryNode(d.pos(i), d.type(i), d.normal(i));
        } else {
            addInternalNode(d.pos(i), d.type(i));
        }
    }
    shape_ = shape() + d.shape();
    return *this;
}

template <class vec_t>
DomainDiscretization<vec_t>& DomainDiscretization<vec_t>::subtract(const DomainDiscretization& d) {
    // Domain 1 is chopped off, and the relevant part of the boundary of `d` is added.
    Range<int> to_add;
    for (int i = 0; i < d.size(); ++i) {
        if (d.type(i) < 0 && shape().contains(d.pos(i))) {
            to_add.push_back(i);
        }
    }
    chopOff(d, to_add);
    for (int i : to_add) {
        addBoundaryNode(d.pos(i), d.type(i), -d.normal(i));
    }
    shape_ = shape() - d.shape();
    return *this;
}

template <class vec_t>
bool DomainDiscretization<vec_t>::contains(const vec_t& point) const {
    assert_msg((shape_->hasContains()), "Domain shape does not have a `contains` method, "
                                        "try using `discreteContains`.");
    return shape_->contains(point);
}

template <class vec_t>
template <typename search_structure_t>
bool DomainDiscretization<vec_t>::discreteContains(
        const vec_t& point, const search_structure_t& search) const {
    assert_msg((search.size() != 0), "Search structure should not be empty.");
    int closest_index = search.query(point, 1).first[0];
    auto closest_point = search.get(closest_index);
    auto n = normals_[closest_index];

    return n.dot(point - closest_point) < 0;
}

template <class vec_t>
template <typename search_structure_t>
bool DomainDiscretization<vec_t>::contains(const vec_t& p, const search_structure_t& search) const {
    if (shape_->hasContains()) return shape_->contains(p);
    else return discreteContains(p, search);
}

template <class vec_t>
template <typename search_structure_t>
void DomainDiscretization<vec_t>::makeDiscreteContainsStructure(search_structure_t& search) const {
    int size = boundary().size();
    Range<vec_t> boundary_points(size);
    for (int i = 0; i < boundary_map_.size(); i++) {
        if (boundary_map_[i] == -1) continue;
        boundary_points[boundary_map_[i]] = positions_[i];
    }
    search.reset(boundary_points);
}

template <class vec_t>
template <typename hdf5_file>
DomainDiscretization<vec_t> DomainDiscretization<vec_t>::load(hdf5_file& file,
                                                              const std::string& name)  {
    file.openGroup(name);
    auto pos = file.readDouble2DArray("pos");
    int n = pos.size();
    UnknownShape<vec_t> shape;
    DomainDiscretization<vec_t> domain(shape);
    domain.positions_.resize(n);
    vec_t low = 1e100;
    vec_t high = -1e100;
    for (int i = 0; i < n; ++i) {
        assert_msg(pos[i].size() == dim, "Node %d has invalid number of coordinates: %s, "
                                         "expected %d.", i, pos[i].size(), dim);
        for (int j = 0; j < dim; ++j) {
            domain.positions_[i][j] = pos[i][j];
            if (pos[i][j] <= low[j]) low[j] = pos[i][j];
            if (pos[i][j] >= high[j]) high[j] = pos[i][j];
        }
    }
    domain.types_ = file.readIntArray("types");
    int bnd_count = (domain.types_ < 0).size();
    domain.boundary_map_ = file.readIntArray("bmap");
    assert_msg(static_cast<int>((domain.boundary_map_ == -1).size()) == n - bnd_count,
               "Number of nodes with boundary index (%d) in bmap is not the same as number "
               "of nodes with negative type (%d).",
               (domain.boundary_map_ == -1).size(), n - bnd_count);
    auto normals = file.readDouble2DArray("normals");
    assert_msg(static_cast<int>(normals.size()) == bnd_count,
            "Number of normals (%d) must be equal to number of boundary nodes (%d).",
            normals.size(), bnd_count);
    domain.normals_.resize(normals.size());
    for (int i = 0; i < bnd_count; ++i) {
        assert_msg(normals[i].size() == dim, "Normal %d has invalid number of coordinates: %s, "
                                             "expected %d.", i, normals[i].size(), dim);
        for (int j = 0; j < dim; ++j) {
            domain.normals_[i][j] = normals[i][j];
        }
    }
    // boundary map consistency check
    for (int i = 0; i < n; ++i) {
        if (domain.types_[i] < 0) {
            assert_msg(0 <= domain.boundary_map_[i] && domain.boundary_map_[i] < bnd_count,
                       "Normals list and boundary map inconsistent. Node %d with type %d "
                       "has boundary index %d, but the list of normals is only %d long.",
                       i, domain.types_[i], domain.boundary_map_[i], bnd_count);
        } else {
            assert_msg(domain.boundary_map_[i] == -1,
                       "Boundary map assigned a normal to an internal node %d with type %d, "
                       "bmap[%d] = %d, but should be -1.",
                       i, domain.types_[i], i, domain.boundary_map_[i]);
        }
    }
    domain.support_.resize(n);
    return domain;
}

template <class vec_t>
DomainDiscretization<vec_t>& DomainDiscretization<vec_t>::translate(const vec_t& a) {
    shape_.reset(new TranslatedShape<vec_t>(*shape_, a));
    for (int i = 0; i < size(); ++i) {
        positions_[i] += a;
    }
    return *this;
}

template <class vec_t>
DomainDiscretization<vec_t>& DomainDiscretization<vec_t>::rotate(
        const Eigen::Matrix<scalar_t, dim, dim>& Q) {
    shape_.reset(new RotatedShape<vec_t>(*shape_, Q));
    for (int i = 0; i < size(); ++i) {
        positions_[i] = Q*positions_[i];
        if (boundary_map_[i] != -1) normals_[boundary_map_[i]] = Q*normals_[boundary_map_[i]];
    }
    return *this;
}

template <class vec_t>
DomainDiscretization<vec_t>& DomainDiscretization<vec_t>::rotate(scalar_t angle) {
    assert_msg(dim == 2, "Angle rotation only available in 2D.");
    scalar_t s = std::sin(angle);
    scalar_t c = std::cos(angle);
    Eigen::Matrix<double, dim, dim> Q; Q << c, -s, s, c;
    return rotate(Q);
}

template <class vec_t>
template <typename compare_t>
Range<int> DomainDiscretization<vec_t>::reorderNodes(const compare_t& cmp) {
    int N = size();
    std::vector<std::pair<vec_t, int>> nodes(N);
    for (int i = 0; i < N; ++i) nodes[i] = {pos(i), i};
    std::sort(nodes.begin(), nodes.end(), cmp);
    indexes_t permutation(N);
    for (int i = 0; i < N; ++i) permutation[i] = nodes[i].second;
    shuffleNodes(permutation);
    return permutation;
}

template <class vec_t>
void DomainDiscretization<vec_t>::shuffleNodes(const indexes_t& permutation) {
    int N = size();
    assert_msg(static_cast<int>(permutation.size()) == N,
            "Permutation size %d must be equal to domain size %d.", permutation.size(), N);
    indexes_t pinv(N, -1);
    for (int i = 0; i < N; ++i) {
        assert_msg(0 <= permutation[i] && permutation[i] < N, "Invalid index %d in permutation.",
                permutation[i]);
        assert_msg(pinv[permutation[i]] == -1, "Duplicate index %d in permutation.",
                permutation[i]);
        pinv[permutation[i]] = i;
    }
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < supportSize(i); ++j) {
            support_[i][j] = pinv[support_[i][j]];
        }
    }
    positions_ = positions_[permutation].asRange();
    types_ = types_[permutation].asRange();
    boundary_map_ = boundary_map_[permutation].asRange();
    support_ = support_[permutation].asRange();
}

}  // namespace mm

#endif  // MEDUSA_BITS_DOMAINS_DOMAINDISCRETIZATION_HPP_
