#ifndef MEDUSA_BITS_APPROXIMATIONS_JACOBISVDWRAPPER_FWD_HPP_
#define MEDUSA_BITS_APPROXIMATIONS_JACOBISVDWRAPPER_FWD_HPP_

/**
 * @file
 * Definition of JacobiSVDWrapper.
 *
 * @example test/approximations/JacobiSVDWrapper_test.cpp
 */

#include <medusa/Config.hpp>
#include <Eigen/Core>
#include <Eigen/SVD>

namespace mm {

/**
 * Extends Eigen's JacobiSVD to compute thin `U` and thin `V` by default.
 * This class exists because is is currently impossible to specify additional
 * parameters to JacobiSVD before calling the compute method.
 * This class satisfies the @ref linsolve-concept.
 * can as such be supplied to approximation engines.
 * @tparam scalar_t Numerical scalar type used in computations.
 * @tparam QRPreconditioner this optional parameter allows to specify the type of QR
 * decomposition that will be used internally for the R-SVD step for non-square matrices.
 * See Eigen docs on JacobiSVD for possible values.
 *
 * Usage example:
 * @snippet approximations/JacobiSVDWrapper_test.cpp JacobiSVDWrapper usage example
 * @ingroup approximations
 */
template <typename scalar_t, int QRPreconditioner = Eigen::ColPivHouseholderQRPreconditioner>
class JacobiSVDWrapper :
        public Eigen::JacobiSVD<Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>,
                                QRPreconditioner> {
  public:
    /**
     * Override compute method by supplying `Eigen::ComputeThinU | Eigen::ComputeThinV`
     * by default.
     */
    void compute(const Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>& M);
};

}  // namespace mm

#endif  // MEDUSA_BITS_APPROXIMATIONS_JACOBISVDWRAPPER_FWD_HPP_
