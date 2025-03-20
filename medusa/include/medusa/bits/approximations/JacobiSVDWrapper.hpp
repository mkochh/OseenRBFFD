#ifndef MEDUSA_BITS_APPROXIMATIONS_JACOBISVDWRAPPER_HPP_
#define MEDUSA_BITS_APPROXIMATIONS_JACOBISVDWRAPPER_HPP_

/**
 * @file
 * Implementation of JacobiSVDWrapper.
 */

#include "JacobiSVDWrapper_fwd.hpp"

namespace mm {

template <typename scalar_t, int QRPreconditioner>
void JacobiSVDWrapper<scalar_t, QRPreconditioner>::compute(
        const Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>& M) {
    Eigen::JacobiSVD<Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>,
                     QRPreconditioner>::compute(
            M, Eigen::ComputeThinU | Eigen::ComputeThinV);
}

}  // namespace mm

#endif  // MEDUSA_BITS_APPROXIMATIONS_JACOBISVDWRAPPER_HPP_
