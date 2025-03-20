/**
 * @file
 * Instantiations of commonly used WLS approximations.
 */

#include <medusa/bits/approximations/WLS.hpp>

#include <medusa/bits/approximations/Gaussian.hpp>
#include <medusa/bits/approximations/Multiquadric.hpp>
#include <medusa/bits/approximations/RBFBasis.hpp>
#include <medusa/bits/approximations/Monomials.hpp>
#include <medusa/bits/approximations/WeightFunction.hpp>
#include <medusa/bits/approximations/ScaleFunction.hpp>
#include <medusa/bits/approximations/JacobiSVDWrapper.hpp>
#include <medusa/bits/types/Vec.hpp>

#include <Eigen/Cholesky>
#include <Eigen/LU>
#include <Eigen/SVD>

namespace mm {
/// Macro for inserting a comma ",".
#define COMMA ,

/// Macro for instantiating template classes with getShape method.
#define INSTANTIATE(classname)  \
    template class classname;   \
    template Eigen::Matrix<double, Eigen::Dynamic, 1> classname::getShape(const Der1<dim>&) const; \
    template Eigen::Matrix<double, Eigen::Dynamic, 1> classname::getShape(const Der2<dim>&) const; \
    template Eigen::Matrix<double, Eigen::Dynamic, 1> classname::getShape(const Lap<dim>&) const; \

INSTANTIATE(WLS<Monomials<Vec2d>>)
INSTANTIATE(WLS<Gaussians<Vec2d>>)

INSTANTIATE(WLS<Monomials<Vec1d> COMMA NoWeight<Vec1d> COMMA ScaleToFarthest>)
INSTANTIATE(WLS<Monomials<Vec2d> COMMA NoWeight<Vec2d> COMMA ScaleToFarthest>)
INSTANTIATE(WLS<Monomials<Vec3d> COMMA NoWeight<Vec3d> COMMA ScaleToFarthest>)

INSTANTIATE(WLS<Monomials<Vec1d> COMMA GaussianWeight<Vec1d> COMMA ScaleToFarthest>)
INSTANTIATE(WLS<Monomials<Vec2d> COMMA GaussianWeight<Vec2d> COMMA ScaleToFarthest>)
INSTANTIATE(WLS<Monomials<Vec2d> COMMA GaussianWeight<Vec2d> COMMA ScaleToClosest>)
INSTANTIATE(WLS<Monomials<Vec3d> COMMA GaussianWeight<Vec3d> COMMA ScaleToClosest>)

INSTANTIATE(WLS<Monomials<Vec2d> COMMA GaussianWeight<Vec2d>>)
INSTANTIATE(WLS<Gaussians<Vec2d> COMMA GaussianWeight<Vec2d>>)
INSTANTIATE(WLS<Gaussians<Vec2d> COMMA GaussianWeight<Vec2d> COMMA ScaleToFarthest>)
INSTANTIATE(WLS<Gaussians<Vec2d> COMMA GaussianWeight<Vec2d> COMMA ScaleToClosest>)
INSTANTIATE(WLS<Gaussians<Vec3d> COMMA GaussianWeight<Vec3d> COMMA ScaleToClosest>)
INSTANTIATE(WLS<Gaussians<Vec3d> COMMA GaussianWeight<Vec3d> COMMA ScaleToFarthest>)
INSTANTIATE(WLS<Gaussians<Vec3d> COMMA NoWeight<Vec3d> COMMA ScaleToFarthest>)

INSTANTIATE(WLS<Gaussians<Vec2d> COMMA NoWeight<Vec2d> COMMA
                    NoScale COMMA Eigen::PartialPivLU<Eigen::MatrixXd>>)
INSTANTIATE(WLS<Gaussians<Vec2d> COMMA NoWeight<Vec2d> COMMA
                    ScaleToFarthest COMMA Eigen::PartialPivLU<Eigen::MatrixXd>>)
INSTANTIATE(WLS<Monomials<Vec2d> COMMA NoWeight<Vec2d> COMMA
                    ScaleToFarthest COMMA Eigen::PartialPivLU<Eigen::MatrixXd>>)

INSTANTIATE(WLS<MQs<Vec2d> COMMA GaussianWeight<Vec2d> COMMA ScaleToClosest>)
INSTANTIATE(WLS<MQs<Vec3d> COMMA GaussianWeight<Vec3d> COMMA ScaleToClosest>)

INSTANTIATE(WLS<Gaussians<Vec2d> COMMA NoWeight<Vec2d> COMMA ScaleToFarthest COMMA
                    Eigen::LLT<Eigen::MatrixXd>>)
INSTANTIATE(WLS<Gaussians<Vec3d> COMMA NoWeight<Vec3d> COMMA ScaleToFarthest COMMA
                    Eigen::LLT<Eigen::MatrixXd>>)
INSTANTIATE(WLS<Gaussians<Vec2d> COMMA NoWeight<Vec2d> COMMA ScaleToFarthest>)

}  // namespace mm
