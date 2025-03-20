/**
 * @file
 * Instantiations of common RBFFD approximations.
 */

#include <medusa/bits/approximations/RBFFD.hpp>

#include <medusa/bits/approximations/Gaussian.hpp>
#include <medusa/bits/approximations/ScaleFunction.hpp>
#include <medusa/bits/approximations/Multiquadric.hpp>
#include <medusa/bits/approximations/InverseMultiquadric.hpp>
#include <medusa/bits/approximations/Polyharmonic.hpp>
#include <medusa/bits/approximations/Monomials.hpp>
#include <medusa/bits/approximations/RBFBasis.hpp>

#include <medusa/bits/types/Vec.hpp>
#include <Eigen/LU>

namespace mm {

/// Macro for inserting a comma ",".
#define COMMA ,
/// Macro for instantiating template classes with getShape method.
#define INSTANTIATE(classname)  \
    template class classname;   \
    template Eigen::Matrix<double , Eigen::Dynamic, 1> classname::getShape(const Der1<dim>&) const;\
    template Eigen::Matrix<double , Eigen::Dynamic, 1> classname::getShape(const Der2<dim>&) const;\
    template Eigen::Matrix<double , Eigen::Dynamic, 1> classname::getShape(const Lap<dim>&) const; \

INSTANTIATE(RBFFD<Gaussian<double> COMMA Vec1d COMMA ScaleToClosest>)
INSTANTIATE(RBFFD<Polyharmonic<double COMMA 3> COMMA Vec1d COMMA ScaleToClosest>)

INSTANTIATE(RBFFD<InverseMultiquadric<double> COMMA Vec2d COMMA ScaleToClosest>)
INSTANTIATE(RBFFD<InverseMultiquadric<double> COMMA Vec3d COMMA ScaleToClosest>)

INSTANTIATE(RBFFD<Multiquadric<double> COMMA Vec2d COMMA ScaleToClosest>)
INSTANTIATE(RBFFD<Multiquadric<double> COMMA Vec3d COMMA ScaleToClosest>)

INSTANTIATE(RBFFD<Gaussian<double> COMMA Vec2d COMMA ScaleToClosest>)
INSTANTIATE(RBFFD<Gaussian<double> COMMA Vec3d COMMA ScaleToClosest>)

INSTANTIATE(RBFFD<Polyharmonic<double COMMA -1> COMMA Vec2d COMMA ScaleToClosest>)
INSTANTIATE(RBFFD<Polyharmonic<double COMMA 3> COMMA Vec2d COMMA ScaleToClosest>)
INSTANTIATE(RBFFD<Polyharmonic<double COMMA 5> COMMA Vec2d COMMA ScaleToClosest>)
INSTANTIATE(RBFFD<Polyharmonic<double COMMA 7> COMMA Vec2d COMMA ScaleToClosest>)

INSTANTIATE(RBFFD<Polyharmonic<double COMMA -1> COMMA Vec3d COMMA ScaleToClosest>)
INSTANTIATE(RBFFD<Polyharmonic<double COMMA 3> COMMA Vec3d COMMA ScaleToClosest>)
INSTANTIATE(RBFFD<Polyharmonic<double COMMA 5> COMMA Vec3d COMMA ScaleToClosest>)
INSTANTIATE(RBFFD<Polyharmonic<double COMMA 7> COMMA Vec3d COMMA ScaleToClosest>)

}  // namespace mm
