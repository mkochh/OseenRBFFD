#ifndef EXAMPLES_CAHNHILLIARD_EQUATION_BIHARMONICOP_HPP_
#define EXAMPLES_CAHNHILLIARD_EQUATION_BIHARMONICOP_HPP_

#include <medusa/Medusa.hpp>

/// Represents the Biharmonic operator (see custom operator example)
template <int dim>
struct Biharmonic : public mm::Operator<Biharmonic<dim>> {
    typedef mm::Vec<double, dim> vec;
    static std::string type_name() { return mm::format("Biharmonic<%d>", dim); }
    static std::string name() { return type_name(); }

    template <int k>
    double applyAt0(const mm::RBFBasis<mm::Polyharmonic<double, k>, vec>&, int index,
                    const std::vector<vec>& support, double scale) const {
        static_assert(k != -1, "If dynamic k is desired it can be obtained from basis.rbf().");
        double r = support[index].norm();
        return k*(k-2)*(dim+k-2)*(dim+k-4)*mm::ipow<k-4>(r) / mm::ipow<4>(scale);
    }

    double applyAt0(const mm::Monomials<vec>& mon, int idx,
                    const std::vector<vec>& q, double s) const {
        double result = 0;
        std::array<int, dim> orders;
        for (int d = 0; d < dim; ++d) { orders[d] = 0; }
        for (int d = 0; d < dim; ++d) {
            orders[d] = 4;
            result += mon.evalOpAt0(idx, mm::Derivative<dim>(orders), q);
            orders[d] = 0;

            for (int d2 = 0; d2 < d; ++d2) {
                orders[d] = 2;
                orders[d2] = 2;
                result += 2*mon.evalOpAt0(idx, mm::Derivative<dim>(orders), q);
                orders[d] = 0;
                orders[d2] = 0;
            }
        }
        return result / mm::ipow<4>(s);
    }
};

#endif  // EXAMPLES_CAHNHILLIARD_EQUATION_BIHARMONICOP_HPP_
