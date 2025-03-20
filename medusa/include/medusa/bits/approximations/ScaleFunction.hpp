#ifndef MEDUSA_BITS_APPROXIMATIONS_SCALEFUNCTION_HPP_
#define MEDUSA_BITS_APPROXIMATIONS_SCALEFUNCTION_HPP_

/**
 * @file
 * Scale function definition and implementation.
 *
 * @example test/approximations/ScaleFunction_test.cpp
 */

#include <medusa/Config.hpp>
#include <medusa/bits/types/Range.hpp>
#include <iosfwd>

namespace mm {

/**
 * Scale function that scales to the closest neighbor. Used as a template parameter
 * in WLS computations of shape functions. This function satisfies the @ref scale-concept.
 * @sa NoScale, ScaleToFarthest
 * @ingroup approximations
 */
class ScaleToClosest {
  public:
    /// Returns local scale of given point and its support.
    template<typename vec_t>
    static inline typename vec_t::scalar_t scale(const vec_t& p,
                                                 const std::vector<vec_t>& support) {
        if (p == support[0]) {
            return (p - support[1]).norm();
        } else {
            return (p - support[0]).norm();
        }
    }
    /// Print information about this scale function.
    inline friend std::ostream& operator<<(std::ostream& os, const ScaleToClosest& /* scale */) {
        return os << "ScaleToClosest";
    }
};

/**
 * Scale function that scales to the farthest neighbor. Used as a template parameter
 * in WLS computations of shape functions. This function satisfies the @ref scale-concept.
 * @sa NoScale, ScaleToClosest
 * @ingroup approximations
 */
class ScaleToFarthest {
  public:
    /// Returns local scale of given point and its support.
    template<typename vec_t>
    static inline typename vec_t::scalar_t scale(const vec_t& p,
                                                 const std::vector<vec_t>& support) {
        return (p - support.back()).norm();
    }
    /// Print information about this scale function.
    inline friend std::ostream& operator<<(std::ostream& os, const ScaleToFarthest& /* scale */) {
        return os << "ScaleToFarthest";
    }
};

/**
 * Scale function that indicates no scaling is performed. Used as a template parameter
 * in WLS computations of shape functions. This function satisfies the @ref scale-concept.
 * @sa ScaleToFarthest, ScaleToClosest
 * @ingroup approximations
 */
class NoScale {
  public:
    template<typename vec_t>
    /// Returns local scale of given point and its support.
    static inline typename vec_t::scalar_t scale(
            const vec_t& /* p */, const std::vector<vec_t>& /* support */) { return 1.0; }
    /// Print information about this scale function.
    inline friend std::ostream& operator<<(std::ostream& os, const NoScale& /* scale */) {
        return os << "NoScale";
    }
};

}  // namespace mm

#endif  // MEDUSA_BITS_APPROXIMATIONS_SCALEFUNCTION_HPP_
