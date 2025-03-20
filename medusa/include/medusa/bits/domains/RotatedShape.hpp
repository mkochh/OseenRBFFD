#ifndef MEDUSA_BITS_DOMAINS_ROTATEDSHAPE_HPP_
#define MEDUSA_BITS_DOMAINS_ROTATEDSHAPE_HPP_

/**
 * @file
 * Implementation of rotated domain shapes.
 */

#include "RotatedShape_fwd.hpp"
#include "DomainDiscretization.hpp"
#include "DomainShape.hpp"
#include <medusa/bits/utils/assert.hpp>

namespace mm {

template <typename vec_t>
RotatedShape<vec_t>::RotatedShape(const DomainShape<vec_t>& sh,
                                  const Eigen::Matrix<scalar_t, dim, dim>& Q) : sh(sh), Q(Q) {
    auto* rsh = dynamic_cast<const RotatedShape<vec_t>*>(&sh);
    if (rsh != nullptr) {  // collapse double rotations
        this->sh = rsh->sh;
        this->Q = this->Q * rsh->Q;
    }
    assert_msg((Q*Q.transpose()).isApprox(Q.Identity(), 1e-12),
            "Matrix Q is not orthogonal, QQ^T - I = %s", (Q*Q.transpose())-Q.Identity());
}

template <typename vec_t>
std::pair<vec_t, vec_t> RotatedShape<vec_t>::bbox() const  {
    auto bb = sh->bbox();
    vec_t m = Q*bb.first;
    vec_t M = m, t;
    for (unsigned i = 1; i < (1 << dim); ++i) {
        for (unsigned j = 0; j < dim; ++j) {
            t[j] = (i & (1u << j)) ? bb.second[j] : bb.first[j];
        }
        t = Q*t;
        for (unsigned j = 0; j < dim; ++j) {
            if (t[j] < m[j]) m[j] = t[j];
            if (t[j] > M[j]) M[j] = t[j];
        }
    }
    return {m, M};
}

/// @cond Doxgen parses this wrong, don't know why...
template <typename vec_t>
DomainDiscretization <vec_t>
RotatedShape<vec_t>::discretizeBoundaryWithStep(scalar_t step, int type) const {
    auto d = sh->discretizeBoundaryWithStep(step, type);
    d.rotate(Q);
    return d;
}

template <typename vec_t>
DomainDiscretization<vec_t>
RotatedShape<vec_t>::discretizeWithStep(scalar_t step, int internal_type, int boundary_type) const {
    auto d = sh->discretizeWithStep(step, internal_type, boundary_type);
    d.rotate(Q);
    return d;
}

template <typename vec_t>
DomainDiscretization<vec_t>
RotatedShape<vec_t>::discretizeWithDensity(const std::function<scalar_t(vec_t)>& dr,
                                           int internal_type, int boundary_type) const {
    auto d = sh->discretizeWithDensity(dr, internal_type, boundary_type);
    d.rotate(Q);
    return d;
}

template <typename vec_t>
DomainDiscretization<vec_t>
RotatedShape<vec_t>::discretizeBoundaryWithDensity(const std::function<scalar_t(vec_t)>& dr,
                                                   int type) const {
    auto d = sh->discretizeBoundaryWithDensity(dr, type);
    d.rotate(Q);
    return d;
}
/// @endcond

}  // namespace mm

#endif  // MEDUSA_BITS_DOMAINS_ROTATEDSHAPE_HPP_
