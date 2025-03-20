#ifndef MEDUSA_BITS_DOMAINS_POLYHEDRONSHAPE_HPP_
#define MEDUSA_BITS_DOMAINS_POLYHEDRONSHAPE_HPP_

#include "PolyhedronShape_fwd.hpp"
#include "DomainDiscretization.hpp"
#include <medusa/bits/spatial_search/KDTreeMutable.hpp>

namespace mm {

template <typename vec_t>
DomainDiscretization<vec_t> PolyhedronShape<vec_t>::discretizeBoundaryWithDensity(
        const std::function<scalar_t(vec_t)>& dr, int type) const {
    #ifdef USE_CGAL6
    // ugly fix to work with cgal 6.x (i don't know how to work with std::optional properly)
    std::optional<Mesh::Property_map<VI, Vector>> opt_vnormals;
    std::optional<Mesh::Property_map<FI, Vector>> opt_fnormals;
    std::optional<Mesh::Property_map<FI, CGAL::Color>> opt_fcolors;
    Mesh::Property_map<VI, Vector> vnormals;
    Mesh::Property_map<FI, Vector> fnormals;
    Mesh::Property_map<FI, CGAL::Color> fcolors;
    opt_vnormals = surface.property_map<VI, Vector>("v:normals");
    assert_msg(opt_vnormals.has_value(), "Vertex normals were not computed.");
    opt_fnormals = surface.property_map<FI, Vector>("f:normals");
    assert_msg(opt_fnormals.has_value(), "Face normals were not computed.");
    opt_fcolors = surface.property_map<FI, CGAL::Color>("f:color");
    vnormals = *std::move(opt_vnormals);
    fnormals = *std::move(opt_fnormals);
    bool cexists = opt_fcolors.has_value();
    fcolors = *std::move(opt_fcolors);
    #endif
    #ifdef USE_CGAL_LESS_THAN_6
    bool vexists, fexists, cexists;
    Mesh::Property_map<VI, Vector> vnormals;
    Mesh::Property_map<FI, Vector> fnormals;
    Mesh::Property_map<FI, CGAL::Color> fcolors;
    std::tie(vnormals, vexists) = surface.property_map<VI, Vector>("v:normals");
    assert_msg(vexists, "Vertex normals were not computed.");
    std::tie(fnormals, fexists) = surface.property_map<FI, Vector>("f:normals");
    assert_msg(fexists, "Face normals were not computed.");
    std::tie(fcolors, cexists) = surface.property_map<FI, CGAL::Color>("f:color");
    #endif
    auto get_face_type = [&](FI face) {
        return (cexists && type == 0) ?
               rgb2type(fcolors[face].r(), fcolors[face].g(), fcolors[face].b()) :
               -1;
    };
    KDTreeMutable<vec_t> tree;
    auto d = mm::DomainDiscretization<vec_t>(*this);
    auto add_point = [&](const vec_t& p, int type, const vec_t& n) {
        if (tree.existsPointInSphere(p, dr(p))) return;
        tree.insert(p);
        d.addBoundaryNode(p, type, n);
    };

    // Discretize vertices.
    for (VI idx : surface.vertices()) {
        int face_type = get_face_type(surface.face(surface.halfedge(idx)));
        auto n = vnormals[idx];
        add_point(getPoint(idx), face_type, {n.x(), n.y(), n.z()});
    }
    // Discretize edges.
    for (Mesh::edge_index idx : surface.edges()) {
        FI face1 = surface.face(idx.halfedge());
        FI face2 = surface.face(surface.opposite(idx.halfedge()));
        Vector normal1 = fnormals[face1];
        Vector normal2 = fnormals[face2];
        vec_t n1 = {normal1.x(), normal1.y(), normal1.z()};
        vec_t n2 = {normal2.x(), normal2.y(), normal2.z()};
        vec_t n = (n1 + n2).normalized();
        int face_type = get_face_type(face1);

        vec_t p1 = getPoint(surface.vertex(idx, 0));
        vec_t p2 = getPoint(surface.vertex(idx, 1));
        for (const vec_t& p : discretization_helpers::discretizeLineWithDensity(p1, p2, dr)) {
            add_point(p, face_type, n);
        }
    }
    // Discretize faces.
    for (FI face : surface.faces()) {
        auto range = surface.vertices_around_face(surface.halfedge(face));
        assert_msg(range.size() == 3, "All faces must be triangles.");
        auto it = range.begin();
        VI p1 = *it;
        ++it;
        VI p2 = *it;
        ++it;
        VI p3 = *it;

        int face_type = get_face_type(face);
        Vector normal = fnormals[face];
        vec_t n = {normal.x(), normal.y(), normal.z()};
        std::vector<vec_t> points_on_surface =
                discretization_helpers::discretizeTriangleWithDensity(
                        getPoint(p1), getPoint(p2), getPoint(p3), n, dr);
        for (const auto& p : points_on_surface) {
            add_point(p, face_type, n);
        }
    }
    return d;
}

}  // namespace mm

#endif  // MEDUSA_BITS_DOMAINS_POLYHEDRONSHAPE_HPP_
