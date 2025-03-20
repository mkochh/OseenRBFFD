#include <medusa/bits/approximations/Monomials_fwd.hpp>
#include <medusa/bits/approximations/WeightFunction_fwd.hpp>
#include <medusa/bits/approximations/ScaleFunction.hpp>
#include <medusa/bits/operators/RaggedShapeStorage.hpp>
#include <medusa/bits/approximations/WLS_fwd.hpp>
#include <medusa/bits/io/HDF.hpp>
#include <medusa/bits/domains/DomainDiscretization_fwd.hpp>
#include <medusa/bits/domains/BoxShape_fwd.hpp>
#include <medusa/bits/domains/FindClosest.hpp>
#include <medusa/bits/operators/UniformShapeStorage.hpp>
#include <medusa/bits/operators/computeShapes.hpp>
#include <medusa/bits/operators/ExplicitOperators.hpp>
#include <medusa/bits/types/ScalarField.hpp>
#include <medusa/bits/utils/numutils.hpp>
#include "Eigen/LU"
#include <medusa/bits/approximations/JacobiSVDWrapper.hpp>

#include "gtest/gtest.h"

using namespace mm;  // NOLINT(*)
using namespace std;  // NOLINT(*)

TEST(End2end, PoissonExplicit) {
    HDF hdf("test.h5", HDF::DESTROY);
    hdf.close();

    BoxShape<Vec2d> box(0.0, 1.0);
    double step = 0.1;
    DomainDiscretization<Vec2d> domain = box.discretizeWithStep(step);
    domain.findSupport(FindClosest(9));

    hdf.atomic().writeDomain("domain", domain);

    WLS<Monomials<Vec2d>, NoWeight<Vec2d>, ScaleToFarthest>
    wls(Monomials<Vec2d>::tensorBasis(2), {});

    auto storage = domain.computeShapes<sh::lap>(wls);
    auto op = storage.explicitOperators();

    ScalarFieldd s = ScalarFieldd::Ones(domain.size());
    ScalarFieldd s2 = ScalarFieldd::Zero(domain.size());
    s[domain.boundary()] = 0.0;

    double T = 0.1;
    double dt = 1e-5;
    int steps = iceil(T/dt);
    auto interior = domain.interior();
    for (int t = 0; t < steps; ++t) {
        for (int i : interior) {
            s2[i] = s[i] + dt * op.lap(s, i);
        }
        s = s2;
    }

    hdf.reopen();
    hdf.writeDoubleArray("sol", s);
    hdf.close();
}
