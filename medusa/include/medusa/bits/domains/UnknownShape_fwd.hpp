#ifndef MEDUSA_BITS_DOMAINS_UNKNOWNSHAPE_FWD_HPP_
#define MEDUSA_BITS_DOMAINS_UNKNOWNSHAPE_FWD_HPP_

#include <medusa/Config.hpp>
#include "DomainShape_fwd.hpp"

/**
 * @file
 * Declaration of UnknownShape class.
 *
 * @example test/domains/UnknownShape_test.cpp
 */

namespace mm {

/**
 * This class represents an unknown domain shape. It can be used for domains when shape is
 * not known or important, like in simple examples:
 * @snippet domains/DomainDiscretization_test.cpp UnknownShape usage example
 *
 * This is the shape used as a placeholder
 * when loading domain discretizations.
 * @snippet domains/DomainDiscretization_test.cpp load usage example
 * @ingroup domains
 */
template <typename vec_t>
class UnknownShape : public DomainShape<vec_t> {
  public:
    using typename DomainShape<vec_t>::scalar_t;
    using DomainShape<vec_t>::dim;

    bool contains(const vec_t& /* point */) const override { return true; }

    bool hasContains() const override { return false; }

    std::pair<vec_t, vec_t> bbox() const override {
        assert_msg(false, "This function is not available for this shape.");
        return {};
    }

    DomainDiscretization<vec_t> discretizeBoundaryWithDensity(
            const std::function<scalar_t(vec_t)>&, int) const override {
        assert_msg(false, "This function is not available for this shape.");
        return DomainDiscretization<vec_t>(*this);
    }

    std::ostream& print(std::ostream& ostream) const override {
        return ostream << "Unknown shape";
    }

    UnknownShape* clone() const override {
        return new UnknownShape<vec_t>(*this);
    }
};

}  // namespace mm

#endif  // MEDUSA_BITS_DOMAINS_UNKNOWNSHAPE_FWD_HPP_
