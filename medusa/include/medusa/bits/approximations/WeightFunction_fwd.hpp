#ifndef MEDUSA_BITS_APPROXIMATIONS_WEIGHTFUNCTION_FWD_HPP_
#define MEDUSA_BITS_APPROXIMATIONS_WEIGHTFUNCTION_FWD_HPP_

/**
 * @file
 * Declarations of weight functions.
 *
 * @example test/approximations/WeightFunction_test.cpp
 */

#include <medusa/Config.hpp>
#include <array>
#include <iosfwd>

namespace mm {

/**
 * Class representing no weight function, i.e.\ a constant 1.
 * Used in WLS computations of shape functions. This class satisfies the @ref weight-concept.
 * Usage example: see WLS.
 * @sa RBFWeight, GaussianWeight, MQWeight, IMQWeight, PHWeight
 * @ingroup approximations
 */
template <class vec_t>
class NoWeight {
  public:
    typedef vec_t vector_t;   ///< Vector type.
    typedef typename vector_t::scalar_t scalar_t;   ///< Scalar type.
    /// Store dimension of the domain.
    enum { /** Dimensionality of the function domain. */ dim = vec_t::dim };

  public:
    /// Evaluate weight function at `point`. Returns `1`.
    scalar_t operator()(const vector_t& /* point */) const { return 1.0; }

    /// Output basic info about given weight function.
    template <typename V>
    friend std::ostream& operator<<(std::ostream& os, const NoWeight<V>& m);
};

/// Output basic info about given weight function.
template <typename V>
std::ostream& operator<<(std::ostream& os, const NoWeight<V>& m) {
    return os << "NoWeight " << m.dim << "D";
}

/**
 * Represents a weight function constructed from a Radial Basis function.
 * @tparam RBFType Type of the RBF used. Must satisfy the @ref rbf-concept.
 * @tparam vec_t Vector type used in calculations.
 *
 * The weight function @f$w@f$ is defined as @f$w(x) = f(\|x\|_2^2)@f$, where @f$f@f$ is the
 * given radial basis function.
 *
 * This class satisfies the @ref weight-concept.
 *
 * Usage example:
 * @snippet approximations/WeightFunction_test.cpp RBF weight usage example
 * @sa NoWeight, GaussianWeight, MQWeight, IMQWeight, PHWeight
 * @ingroup approximations
 */
template <class RBFType, class vec_t>
class RBFWeight {
    static_assert(std::is_same<typename vec_t::scalar_t, typename RBFType::scalar_t>::value,
                  "Basis and underlying RBF must have the same scalar type.");
  public:
    typedef RBFType rbf_t;   ///< Radial basis function type.
    typedef vec_t vector_t;   ///< Vector type.
    typedef typename vector_t::scalar_t scalar_t;   ///< Scalar type.
    /// Store dimension of the domain.
    enum { /** Dimensionality of the function domain. */ dim = vec_t::dim };

  private:
    rbf_t rbf_;  ///< RBF function.

  public:
    /// Construct a weight with default construction of RBF.
    RBFWeight() : rbf_() {}
    /// Construct a weight from given RBF.
    RBFWeight(const rbf_t& rbf) : rbf_(rbf) {}
    /// Perfect forwarding constructor. Construct RBFWeight as if you were constructing given RBF.
    template <typename ...Args>
    RBFWeight(Args&&... args) : rbf_(std::forward<Args>(args)...) {}

    /// Returns underlying RBF object.
    const rbf_t& rbf() const { return rbf_; }

    /// Returns modifiable underlying RBF object.
    rbf_t& rbf() { return rbf_; }

    /// Evaluate weight function at `point`.
    scalar_t operator()(const vector_t& point) const { return rbf_(point.squaredNorm()); }

    /// Output basic info about given weight function.
    template <typename V, typename R>
    friend std::ostream& operator<<(std::ostream& os, const RBFWeight<V, R>& m);
};

// Convenience typedefs

template <typename S> class Gaussian;
/// RBF weight function using Gaussian RBF. Defined for convenience.
template <typename V> using GaussianWeight = RBFWeight<Gaussian<typename V::scalar_t>, V>;

template <typename S> class Multiquadric;
/// RBF weight function using Multiquadric RBF. Defined for convenience.
template <typename V> using MQWeight = RBFWeight<Multiquadric<typename V::scalar_t>, V>;

template <typename S> class InverseMultiquadric;
/// RBF weight function using InverseMultiquadric RBF. Defined for convenience.
template <typename V> using IMQWeight = RBFWeight<InverseMultiquadric<typename V::scalar_t>, V>;

template <typename S, int k> class Polyharmonic;
/// RBF weight function using Polyharmonic RBF. Defined for convenience.
template <typename V> using PHWeight = RBFWeight<Polyharmonic<typename V::scalar_t, -1>, V>;

}  // namespace mm

#endif  // MEDUSA_BITS_APPROXIMATIONS_WEIGHTFUNCTION_FWD_HPP_
