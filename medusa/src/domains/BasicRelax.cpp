#include <medusa/bits/domains/BasicRelax.hpp>

/**
 * @file
 * Instantiation of class for node relaxation.
 */

namespace mm {

BasicRelax& BasicRelax::onlyNodes(Range<int> nodes) {
    nodes_ = nodes; return *this;
}

BasicRelax& BasicRelax::initialHeat(double in) {
    initial_heat = in; return *this;
}

BasicRelax& BasicRelax::finalHeat(double in) {
    final_heat = in; return *this;
}

BasicRelax& BasicRelax::projectionType(ProjectionType in) {
    projection_type = in; return *this;
}

BasicRelax& BasicRelax::boundaryProjectionThreshold(double in) {
    assert_msg(0.0 <= in && in < 2.0, "Weird projection threshold, it is equal to %.6f, "
                                      "but expected in range [0, 1].", in);
    boundary_projection_threshold = in; return *this;
}


BasicRelax& BasicRelax::numNeighbours(int neighbours) {
    num_neighbours = neighbours; return *this;
}

BasicRelax& BasicRelax::iterations(int iterations) {
    num_iterations = iterations; return *this;
}

BasicRelax& BasicRelax::potentialOrder(int order) {
    potential_order = order; return *this;
}

BasicRelax& BasicRelax::rebuildTreeAfter(int iterations) {
    rebuild_tree_after = iterations; return *this;
}

}  // namespace mm
