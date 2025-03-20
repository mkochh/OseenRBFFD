#include <medusa/Medusa.hpp>
#include <medusa/bits/domains/GeneralFill.hpp>

using namespace mm;  // NOLINT
using namespace std;  // NOLINT
using namespace Eigen;  // NOLINT

/// Definition of custom stencil selection algorithms.
/// http://e6.ijs.si/medusa/wiki/index.php/Customization

/// A bad (inefficient and not robust) visibility based stencil finder.
struct VisibilityStencil {
    int num_closest;
    explicit VisibilityStencil(int num_closest) : num_closest(num_closest) {}

    /// Check for straight-line visibility by testing intermediate points with some density.
    template <typename vec_t>
    bool visible(const vec_t& p, const vec_t& q, const DomainDiscretization<vec_t>& domain) const {
        int density = 100;
        vec_t increment = (q-p) / density;
        for (int i = 1; i < density; ++i) {
            if (!domain.contains(p + i*increment)) {
                return false;
            }
        }
        return true;
    }

    /// Method responsible for stencil construction.
    template <typename vec_t>
    void operator()(DomainDiscretization<vec_t>& domain) const {
        int N = domain.size();
        KDTree<vec_t> tree(domain.positions());
        for (int i = 0; i < N; ++i) {
            Range<int> indices = tree.query(domain.pos(i), num_closest).first;
            for (int j : indices) {
                if (visible(domain.pos(i), domain.pos(j), domain)) {
                    domain.support(i).push_back(j);
                }
            }
        }
    }
};


int main() {
    PolygonShape<Vec2d> poly({{0.5, 0.0}, {1.5, 0.0}, {1.5, 0.6}, {1.1, 0.6}, {1.0, 0.2},
                              {0.9, 0.6}, {0.5, 0.6}});
    GeneralFill<Vec2d> fill; fill.seed(0);
    DomainDiscretization<Vec2d> domain = poly.discretizeWithDensity(
            [](const Vec2d&) { return 0.05; }, fill);
    DomainDiscretization<Vec2d> domain_copy = domain;

    int n = 25;
    domain.findSupport(VisibilityStencil(n));
    domain_copy.findSupport(FindClosest(n));

    std::ofstream out_file("custom_stencils_visibility_data.m");
    out_file << "positions = " << domain.positions() << ";" << std::endl;
    out_file << "stencils_closest = " << domain_copy.supports() << ";" << std::endl;
    out_file << "stencils_visibility = " << pad(domain.supports(), -1) << ";" << std::endl;

    return 0;
}

