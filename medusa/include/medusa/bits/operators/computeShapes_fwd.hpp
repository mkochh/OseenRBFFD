#ifndef MEDUSA_BITS_OPERATORS_COMPUTESHAPES_FWD_HPP_
#define MEDUSA_BITS_OPERATORS_COMPUTESHAPES_FWD_HPP_

/**
 * @file
 * Declarations of shape function computation utilities.
 *
 * @example test/operators/computeShapes_test.cpp
 */

#include <medusa/Config.hpp>
#include <medusa/bits/operators/shape_flags.hpp>
#include <iosfwd>

namespace mm {

template <typename vec_t>
class DomainDiscretization;

/**
 * Computes shape functions (stencil weights) for given nodes using support domains from `domain`
 * and approximations provided by `approx`. The computed shapes are stored in `storage`.
 *
 * @tparam approx_t Type of the approximation engine to use.
 * Must satisfy the @ref approx-concept.
 * @tparam shape_storage_t Type of the shape storage. Must satisfy the @ref ss-concept.
 * @tparam OpFamilies A list of operator families for which to compute the shapes. The families must
 * satisfy the @ref operator-family-concept. The family types must be unique.
 * @param[in] domain Domain discretization for which to compute the shape functions.
 * @param[in] approx Approximation engine specifying the shape computation procedure.
 * @param[in] indexes A set of indexes for which to compute the shape functions. Common
 * values are @ref DomainDiscretization::all "domain.all()",
 * @ref DomainDiscretization::interior "domain.interior()",
 * @ref DomainDiscretization::boundary "domain.boundary()".
 * @param[in] operators Object representing operator families.
 * @param[out] storage An object for storing the shape functions. The object is filled
 * during computation, but must be appropriately reserved before computing shapes, by using
 * @ref ShapeStorage::resize.
 *
 * @note This function supports parallel execution if OpenMP is enabled in the compiler, e.g.\
 * for GCC: `-fopenmp`, for ICC: `-openmp`.
 *
 * Usage example:
 * @snippet operators/computeShapes_test.cpp computeShapes usage example
 *
 * Usage with (dummy) custom operator:
 * @snippet operators/computeShapes_test.cpp custom zero op
 * @snippet operators/computeShapes_test.cpp custom op usage
 *
 * @ingroup operators
 */
template <typename approx_t, typename shape_storage_t, typename ...OpFamilies>
void computeShapes(const DomainDiscretization<typename approx_t::vector_t>& domain,
                   approx_t approx, const indexes_t& indexes,
                   const std::tuple<OpFamilies...>& operators, shape_storage_t* storage);

}  // namespace mm

#endif  // MEDUSA_BITS_OPERATORS_COMPUTESHAPES_FWD_HPP_
