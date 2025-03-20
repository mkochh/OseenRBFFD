#include <medusa/Medusa_fwd.hpp>

/// Basic medusa example showcasing what to do, when points "leak" out of the domain
/// https://e6.ijs.si/medusa/wiki/index.php/Determining_the_interior_of_the_domain_by_oversampling_the_boundary

using namespace mm;  // NOLINT

int main() {
    // Points to be used in NURBS creation.
    Range<Vec2d> pts{{210, -132.5}, {205, -115}, {125, -35}, {85, -61}, {85, -61}, {85, -61},
                     {80, -65}, {75, -58}, {75, -58}, {65, 45}, {25, 25}, {-15, -16}, {-15, -16},
                     {-15, -16}, {-35, -28}, {-40, -32.375}, {-40, -32.375}, {-43, -35}, {-27, -40},
                     {-5, -40}, {-5, -40}, {50, -38}, {35, -65}, {15, -75}, {15, -75}, {-15, -95},
                     {-20, -140}, {45, -146}, {45, -146}, {95, -150}, {215, -150}, {210, -132.5}};

    // Calculate NURBS patches' control points, weights, knot vector and order.
    int p = 3;
    Range<double> weights{1, 1, 1, 1};
    Range<double> knot_vector{0, 0, 0, 0, 1, 1, 1, 1};
    Range<NURBSPatch<Vec2d, Vec1d>> patches;

    for (int i = 0; i < 8; i++) {
        Range<Vec2d> control_points;
        for (int j = 0; j < 4; j++) {
            control_points.push_back(pts[4 * i + j]);
        }

        NURBSPatch<Vec2d, Vec1d> patch(control_points, weights, {knot_vector}, {p});
        patches.push_back(patch);
    }

    // Create NURBS shape from a range of patches.
    NURBSShape<Vec2d, Vec1d> shape(patches);
    shape.seed(1);

    // Fill the domain - we get point leakage.
    double h = 5;
    DomainDiscretization<Vec2d> domain_leak = shape.discretizeBoundaryWithStep(h);

    GeneralFill<Vec2d> gf1;
    gf1.seed(1).maxPoints(100000);
    gf1(domain_leak, h);

    // Oversample the domain.
    DomainDiscretization<Vec2d> oversampled_domain = shape.discretizeBoundaryWithStep(h / 5);

    // Construct the contains function.
    KDTree<Vec2d> contains_tree;
    oversampled_domain.makeDiscreteContainsStructure(contains_tree);
    auto contains_function = [&] (const Vec2d p) {
        return oversampled_domain.discreteContains(p, contains_tree);
    };

    // Fill the boundary normally.
    DomainDiscretization<Vec2d> domain = shape.discretizeBoundaryWithStep(h);

    // Fill the interior with the contains function.
    KDTreeMutable<Vec2d> tree;
    GeneralFill<Vec2d> gf;
    gf.seed(1).maxPoints(100000);
    gf(domain, h, tree, contains_function);


    // Write the solution into file.
    std::ofstream out_file("inside_by_oversampling_data.m");
    out_file << "positions_leaking = " << domain_leak.positions() << ";" << std::endl;
    out_file << "positions = " << domain.positions() << ";" << std::endl;
    out_file.close();

    return 0;
}
