#include <medusa/bits/interpolants/PUApproximant.hpp>
#include <medusa/bits/approximations/RBFFD.hpp>
#include <medusa/bits/approximations/Polyharmonic.hpp>

/**
 * @file
 * Instantiation of class for partition of unity approximant.
 */

namespace mm {

template class PUApproximant<Vec1d>;
template class PUApproximant<Vec2d>;
template class PUApproximant<Vec3d>;

/// @cond

// 1D
template Eigen::VectorXd PUApproximant<Vec1d>::evaluate(
    const DomainDiscretization<Vec1d>& domain, const Eigen::VectorXd& values,
    const Range<Vec1d>& query_points, double radius_factor,
    const RBFFD<Polyharmonic<double>, Vec1d, ScaleToClosest>& engine);

template Eigen::VectorXd PUApproximant<Vec1d>::evaluate(
    const DomainDiscretization<Vec1d>& domain, const Eigen::VectorXd& values,
    const KDTree<Vec1d>& query_points_tree, double radius_factor,
    const RBFFD<Polyharmonic<double>, Vec1d, ScaleToClosest>& engine);

// 2D.
template Eigen::VectorXd PUApproximant<Vec2d>::evaluate(
    const DomainDiscretization<Vec2d>& domain, const Eigen::VectorXd& values,
    const Range<Vec2d>& query_points, double radius_factor,
    const RBFFD<Polyharmonic<double>, Vec2d, ScaleToClosest>& engine);

template Eigen::VectorXd PUApproximant<Vec2d>::evaluate(
    const DomainDiscretization<Vec2d>& domain, const Eigen::VectorXd& values,
    const KDTree<Vec2d>& query_points_tree, double radius_factor,
    const RBFFD<Polyharmonic<double>, Vec2d, ScaleToClosest>& engine);

// 3D.
template Eigen::VectorXd PUApproximant<Vec3d>::evaluate(
    const DomainDiscretization<Vec3d>& domain, const Eigen::VectorXd& values,
    const Range<Vec3d>& query_points, double radius_factor,
    const RBFFD<Polyharmonic<double>, Vec3d, ScaleToClosest>& engine);

template Eigen::VectorXd PUApproximant<Vec3d>::evaluate(
    const DomainDiscretization<Vec3d>& domain, const Eigen::VectorXd& values,
    const KDTree<Vec3d>& query_points_tree, double radius_factor,
    const RBFFD<Polyharmonic<double>, Vec3d, ScaleToClosest>& engine);

/// @endcond

}  // namespace mm
