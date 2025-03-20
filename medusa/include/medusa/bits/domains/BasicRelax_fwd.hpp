#ifndef MEDUSA_BITS_DOMAINS_BASICRELAX_FWD_HPP_
#define MEDUSA_BITS_DOMAINS_BASICRELAX_FWD_HPP_

#include <medusa/Config.hpp>
#include <medusa/bits/types/Range_fwd.hpp>

/**
 * @file
 * Declaration of class for domain relaxation.
 *
 * @example test/domains/BasicRelax_test.cpp
 */

namespace mm {


/**
 * Redistributes nodes towards more uniform distribution by minimizing potential between nodes. The
 * engine first finds `num_neighbours` closest neighbours. Then it translates each point according
 * to the repelling "force" induced by its neighbours, where it also takes into account target
 * potential. The force is determined by summing the gradients of potentials in surrounding nodes,
 * where amplitude is determined based on distance to the closest neighbour and supplied target
 * density function. The algorithm uses also a primitive Simulated Annealing, i.e.\ the relax
 * movement magnitude is multiplied with an annealing factor that is linearly dropping from
 * initial_heat to final_heat.
 *
 * The Engine supports two call types: with supplied distribution function relax(func), where it
 * tries to satisfy the user supplied nodal density function. This can be achieved only when there
 * is the total number of domain nodes the same as integral of density function over the domain. If
 * there is too much nodes a volatile relax might occur. If there is not enough nodes the relax
 * might become lazy. The best use of this mode is in combination with fillDistribution Engines,
 * check test for examples.
 *
 * Without distribution, where nodes always move towards less populated area regardless anything.
 * The relax magnitude is simply determined from annealing factor and distance to the closest node.
 * A simple and stable approach, however, note that this relax always converges towards uniformly
 * distributed nodes.
 *
 * Usage example:
 * @snippet domains/BasicRelax_test.cpp BasicRelax usage example
 * @ingroup domains
 */
class BasicRelax {
  public:
    /// Indicating type of projection used when a relax node goes out of the domain.
    enum ProjectionType {
        /**
         * Escaped nodes are frozen outside the domain till the end of relax and
         * then removed. This is good option for violent and unstable relaxations.
         */
        DO_NOT_PROJECT,
        /**
         * Project on boundary in relax movement direction. In general produces better
         * distribution on boundaries, but the approximated normals can be bad.
         */
        PROJECT_IN_DIRECTION,
        /**
         * Project between two closest boundary nodes – this one might result in divergent
         * behaviour, when both closest nodes are on one side, leading in huge gaps on the
         * boundary. However, when there are enough boundary nodes in the first place, this
         * works very well.
         */
        PROJECT_BETWEEN_CLOSEST
    };

  private:
    int num_neighbours = 1;  ///< Number of nodes to consider when calculating the potential.
    int num_iterations = 50;  ///< Number of iterations performed.
    double initial_heat = 1;  ///< Initial heat, usually between 0 and 5.
    double final_heat = 0;  ///< Heat at the end of the relax, usually around 0.
    int potential_order = 2;  ///< Order of repulsing potential.
    int rebuild_tree_after = 1;  ///< How often engine rebuild search tree, 1 is perfect but slow.
    Range<int> nodes_;  ///< List of nodes to process.
    ProjectionType projection_type = DO_NOT_PROJECT;  ///< On boundary projection method.
    double boundary_projection_threshold = 0.75;  ///< Threshold for projecting nodes on boundary.

  public:
    BasicRelax() = default;

    /// Move only given nodes.
    BasicRelax& onlyNodes(Range<int> nodes);

    /**
     * Sets initial heat.
     * @param in Initial heat, usually between 0 and 5, higher heat
     * means more volatile movement. High initial heat may cause divergence and erratic behaviour.
     * Setting too small initial heat might results in lazy relaxation.
     */
    BasicRelax& initialHeat(double in);

    /// Sets final heat.
    BasicRelax& finalHeat(double in);

    /// Sets num neighbours.
    BasicRelax& numNeighbours(int neighbours);

    /// Sets order of repulsing potential
    BasicRelax& potentialOrder(int order);

    /// Sets number of iterations.
    BasicRelax& iterations(int iterations);

    /**
     * Sets rebuild tree frequency. Ir rebuild's tree every `in` iterations.
     * `in = 1` is perfect but slow. Using higher values results in better
     * performance at the cost of accuracy.
     */
    BasicRelax& rebuildTreeAfter(int iterations);

    /// Determines how to handle nodes that escape during relaxation.
    BasicRelax& projectionType(ProjectionType in);

    /**
     * Sets threshold for adding nodes on boundary, i.e.\ if node @f$d_1@f$ and @f$d_2@f$
     * are distances to closest boundary nodes, do not add node on boundary if
     * @f$d_1/d_2 < boundary_projection_threshold @f$.
     * If threshold is 0 all nodes are added, if it is greater than 1 no nodes are added.
     */
    BasicRelax& boundaryProjectionThreshold(double in);

    /**
     * Runs the relax on the selected domain with constant distribution equals
     * to domain characteristic distance
     * @param domain domain to process
     */
    template<class domain_t>
    void operator()(domain_t& domain) const {
        typedef typename domain_t::vector_t vec_t;
        operator()(domain, [](const vec_t& /* p */) { return -1.0; });
    }

    /**
     * Runs the relax on the selected domain with constant density
     * @param domain domain to process
     * @param r constant density
    */
    template<class domain_t>
    void operator()(domain_t& domain, double r) const {
        typedef typename domain_t::vector_t vec_t;
        operator()(domain, [r](const vec_t& /* p */) { return r; });
    }

    /// Runs the procedure of a given domain.
    template<class domain_t, class radius_func_type>
    void operator()(domain_t& domain, const radius_func_type& r_func) const;
};  // class BasicRelax

}  // namespace mm

#endif  // MEDUSA_BITS_DOMAINS_BASICRELAX_FWD_HPP_
