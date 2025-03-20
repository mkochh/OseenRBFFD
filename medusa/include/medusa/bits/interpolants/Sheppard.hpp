#ifndef MEDUSA_BITS_INTERPOLANTS_SHEPPARD_HPP_
#define MEDUSA_BITS_INTERPOLANTS_SHEPPARD_HPP_

/**
 * @file
 * Implementation of SheppardInterpolant class.
 */

#include <medusa/bits/utils/assert.hpp>
#include <medusa/bits/utils/numutils.hpp>

#include "Sheppard_fwd.hpp"

namespace mm {

template <class vec_t, class value_t>
value_t SheppardInterpolant<vec_t, value_t>::operator()(const vec_t& point, int num_closest,
                                                        int power, double reg,
                                                        double conf_dist) const {
    int N = tree.size();
    int M = values.size();
    assert_msg(N == M && N >= 0,
               "Missing or wrong number of positions/values. Did you forget to call setValues() or "
               "setPositions()?");
    assert_msg(num_closest >= 0, "Number of closest nodes should be greater than 0, got %d.",
               num_closest);
    assert_msg(reg >= 0, "Regularization parameter should be greater than 0, got %e.", reg);
    assert_msg(power >= 0, "Power should be greater than 0, got %e.", power);
    assert_msg(num_closest <= N, "You requested %d closest points, but there are only %d in total.",
               num_closest, N);

    Range<value_t> d2;  // Squares of distances.
    Range<int> idx;     // Closest nodes indexes.
    std::tie(idx, d2) = tree.query(point, num_closest);

    int n = idx.size();
    scalar_t s = 0.0, r = 0.0;
    Range<scalar_t> inv_dist(n);

    if (d2[0] < conf_dist * conf_dist) {
        // Point is close to the tree point, thus return value in node.
        return values[idx[0]];
    }

    for (int i = 0; i < n; ++i) {
        inv_dist[i] = 1.0 / (sqrt(ipow(d2[i], power)) + reg);
        s += inv_dist[i];
    }
    for (int i = 0; i < n; ++i) {
        inv_dist[i] /= s;
    }
    for (int i = 0; i < n; ++i) {
        r += inv_dist[i] * values[idx[i]];
    }

    return r;
}
}  // namespace mm

#endif  // MEDUSA_BITS_INTERPOLANTS_SHEPPARD_HPP_
