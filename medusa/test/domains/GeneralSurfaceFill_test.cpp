#include <medusa/bits/domains/GeneralSurfaceFill.hpp>
#include <medusa/bits/types/Vec.hpp>
#include <medusa/bits/types/Range_fwd.hpp>
#include <medusa/bits/domains/BoxShape.hpp>
#include <medusa/bits/domains/DomainDiscretization_fwd.hpp>
#include <medusa/bits/spatial_search/KDTree.hpp>
#include <string>
#include <fstream>

#include "gtest/gtest.h"

namespace mm {

TEST(DomainEngines, GeneralSurfaceFill2DConstantFill) {
    // Define surface.
    auto example_r = [](Vec1d t) {
        t(0) = -t(0);
        double r = pow(abs(cos(1.5 * t(0))), sin(3 * t(0)));
        return Vec2d(r * cos(t(0)), r * sin(t(0)));
    };

    auto der_example_r = [](Vec1d t) {
        t(0) = -t(0);
        double r = pow(abs(cos(1.5 * t(0))), sin(3 * t(0)));
        double der_r = (-1.5 * pow(abs(cos(1.5 * t(0))),
                sin(3 * t(0))) * sin(3 * t(0)) * sin(1.5 * t(0))  +
                        3 * pow(abs(cos(1.5 * t(0))),
                                sin(3 * t(0))) * cos(3 * t(0)) * cos(1.5 * t(0))
                                * log(abs(cos(1.5 * t(0))))) / cos(1.5 * t(0));

        Vec2d jm;
        jm.col(0) << -(der_r * cos(t(0)) - r * sin(t(0))),
                        -(der_r * sin(t(0)) + r * cos(t(0)));

        return jm;
    };

    double h = 0.02;

    UnknownShape<Vec2d> shape;
    DomainDiscretization<Vec2d> domain(shape);

    // Fill domain.
    GeneralSurfaceFill<Vec2d, Vec1d> gsf; gsf.seed(0);  // deterministic
    BoxShape<Vec1d> bs(Vec1d{0.0}, Vec1d{2 * PI});
    domain.fill(gsf, bs, example_r, der_example_r, h);

    // Calculate distances.
    KDTree<Vec2d> tree(domain.positions());
    Range<double> distances;
    double avg_dist = 0;
    for (const auto& pt : domain.positions()) {
        distances.push_back(sqrt(tree.query(pt, 2).second[1]));
        avg_dist += distances[distances.size() - 1];
    }
    avg_dist /= domain.positions().size();

    auto err_avg = std::abs(avg_dist - h) / h;
    EXPECT_LT(err_avg, 0.05);

    auto err_min = std::abs(*std::min_element(distances.begin() + 1, distances.end() - 1) - h) / h;
    EXPECT_LT(err_min, 0.35);
}

TEST(DomainEngines, GeneralSurfaceFill2DFunctionFill) {
    // Define surface.
    auto example_r = [](Vec1d t) {
        t(0) = -t(0);
        double r = pow(abs(cos(1.5 * t(0))), sin(3 * t(0)));
        return Vec2d(r * cos(t(0)), r * sin(t(0)));
    };

    auto der_example_r = [](Vec1d t) {
        t(0) = -t(0);
        double r = pow(abs(cos(1.5 * t(0))), sin(3 * t(0)));
        double der_r = (-1.5 * pow(abs(cos(1.5 * t(0))),
                sin(3 * t(0))) * sin(3 * t(0)) * sin(1.5 * t(0))  +
                        3 * pow(abs(cos(1.5 * t(0))),
                                sin(3 * t(0))) * cos(3 * t(0)) * cos(1.5 * t(0))
                                * log(abs(cos(1.5 * t(0))))) / cos(1.5 * t(0));

        Eigen::Matrix<double, 2, 1> jm;
        jm.col(0) << -(der_r * cos(t(0)) - r * sin(t(0))), -(der_r * sin(t(0)) + r * cos(t(0)));

        return jm;
    };

    auto gradient_h = [](Vec2d p){
        double h_0 = 0.005;
        double h_m = 0.03 - h_0;

        return  0.5 * h_m * (p(0) + p(1) + 3.0) + h_0;
    };

    UnknownShape<Vec2d> shape;
    DomainDiscretization<Vec2d> domain(shape);

    // Fill domain.
    GeneralSurfaceFill<Vec2d, Vec1d> gsf; gsf.seed(0);  // deterministic
    BoxShape<Vec1d> bs(Vec1d{0.0}, Vec1d{2 * PI});
    DomainDiscretization<Vec1d> param_d(bs);
    domain.fill(gsf, param_d, example_r, der_example_r, gradient_h);

    // Calculate distances.
    KDTree<Vec2d> tree(domain.positions());
    auto positions = domain.positions();
    for (int i = 1; i < positions.size() - 1 ; i++) {
        auto pt = positions[i];
        double dist = sqrt(tree.query(pt, 2).second[1]);
        double h = gradient_h(pt);
        double err = std::abs(dist - h) / h;
        EXPECT_LT(err, 0.35);
    }
}

TEST(DomainEngines, GeneralSurfaceFill3DConstantFill) {
    auto torus_r = [](Vec2d t){
        double a = 10.0, b = 25.0;
        return Vec3d((a * cos(t(1)) + b) * cos(t(0)),
                (a * cos(t(1)) + b) * sin(t(0)), a * sin(t(1)));
    };

    auto torus_jacobian_matrix = [](Vec2d t){
        double a = 10.0, b = 25.0;
        Eigen::Matrix<double, 3, 2> jm;
        jm.col(0) << - (a * cos(t(1)) + b) * sin(t(0)), (a * cos(t(1)) + b) * cos(t(0)), 0;
        jm.col(1) << - a * sin(t(1)) * cos(t(0)), - a * sin(t(1)) * sin(t(0)), a * cos(t(1));
        return jm;
    };

    double h = 1.0;

    UnknownShape<Vec3d> shape;
    DomainDiscretization<Vec3d> domain(shape);

    // Fill domain.
    GeneralSurfaceFill<Vec3d, Vec2d> gsf; gsf.seed(0);  // deterministic
    BoxShape<Vec2d> bs(Vec2d{0.0, 0.0}, Vec2d{2 * PI, 2 * PI});
    DomainDiscretization<Vec2d> param_d(bs);
    domain.fill(gsf, param_d, torus_r, torus_jacobian_matrix, h);

    // Calculate distances.
    KDTree<Vec3d> tree(domain.positions());
    Range<double> distances;
    double avg_dist = 0;
    for (const auto& pt : domain.positions()) {
        auto d = tree.query(pt, 5).second;
        double val = 0.0;
        for (int i = 1; i < d.size(); i++) val += sqrt(d[i]);
        val /= (d.size() - 1);

        distances.push_back(val);
        avg_dist += distances[distances.size() - 1];
    }
    avg_dist /= domain.positions().size();

    auto err_avg = std::abs(avg_dist - h) / h;
    EXPECT_LT(err_avg, 0.1);

    auto err_min = std::abs(*std::min_element(distances.begin() + 1, distances.end() - 1) - h) / h;
    EXPECT_LT(err_min, 0.2);
}

TEST(DomainEngines, GeneralSurfaceFill3DFunctionFill) {
    /// [GeneralSurfaceFill 3d function fill usage example]
    // Define parametrization function.
    auto torus_r = [](Vec2d t){
        double a = 10.0, b = 25.0;
        return Vec3d((a * cos(t(1)) + b) * cos(t(0)),
                     (a * cos(t(1)) + b) * sin(t(0)), a * sin(t(1)));
    };

    // Define Jacobian of parametrization function.
    auto torus_jacobian = [](Vec2d t){
        double a = 10.0, b = 25.0;
        Eigen::Matrix<double, 3, 2> jm;
        jm.col(0) << - (a * cos(t(1)) + b) * sin(t(0)), (a * cos(t(1)) + b) * cos(t(0)), 0;
        jm.col(1) << - a * sin(t(1)) * cos(t(0)), - a * sin(t(1)) * sin(t(0)), a * cos(t(1));
        return jm;
    };

    // Define density function.
    auto gradient_h = [](Vec3d p){
        double h_0 = 1.0;
        return h_0 + (p[0] + p[1] + p[2]) / 500.0;
    };

    // Define domain.
    UnknownShape<Vec3d> shape;
    DomainDiscretization<Vec3d> domain(shape);

    // Fill domain.
    GeneralSurfaceFill<Vec3d, Vec2d> gsf; gsf.seed(0);  // deterministic
    BoxShape<Vec2d> param_domain_shape(Vec2d{0.0, 0.0}, Vec2d{2 * PI, 2 * PI});
    DomainDiscretization<Vec2d> param_domain(param_domain_shape);
    domain.fill(gsf, param_domain, torus_r, torus_jacobian, gradient_h);
    /// [GeneralSurfaceFill 3d function fill usage example]

    // Calculate distances.
    KDTree<Vec3d> tree(domain.positions());
    auto positions = domain.positions();
    for (int i = 1; i < positions.size() - 1 ; i++) {
        auto pt = positions[i];
        auto d = tree.query(pt, 7).second;
        double val = 0.0;
        for (int i = 1; i < d.size(); i++) val += sqrt(d[i]);
        val /= (d.size() - 1);

        double h = gradient_h(pt);
        EXPECT_GE(val, 0.75 * h);
    }
}

TEST(DomainEngines, GeneralSurfaceFill4DConstantFill) {
    auto three_sphere_r = [](Vec3d t){
        double R = 10.0;
        return Vec<double, 4>(R * cos(t(0)), R * sin(t(0)) * cos(t(1)),
                R * sin(t(0)) * sin(t(1)) * cos(t(2)),
                R * sin(t(0)) * sin(t(1)) * sin(t(2)));
    };

    auto three_sphere_jacobian_matrix = [](Vec3d t){
        double R = 10.0;
        Eigen::Matrix<double, 4, 3> jm;
        jm.col(0) << - R * sin(t(0)), R * cos(t(0)) * cos(t(1)),
                        R * cos(t(0)) * sin(t(1)) * cos(t(2)),
                        R * cos(t(0)) * sin(t(1)) * sin(t(2));
        jm.col(1) << 0, - R * sin(t(0)) * sin(t(1)),
                        R * sin(t(0)) * cos(t(1)) * cos(t(2)),
                        R * sin(t(0)) * cos(t(1)) * sin(t(2));
        jm.col(2) << 0, 0, - R * sin(t(0)) * sin(t(1)) * sin(t(2)),
                        R * sin(t(0)) * sin(t(1)) * cos(t(2));
        return jm;
    };

    double h = 10.0;

    // Define domain.
    UnknownShape<Vec<double, 4> > shape;
    DomainDiscretization<Vec<double, 4> > domain(shape);

    // Fill domain.
    GeneralSurfaceFill<Vec<double, 4>, Vec3d> gsf; gsf.seed(0);  // deterministic
    BoxShape<Vec3d> bs(Vec3d{0.0, 0.0, 0.0}, Vec3d{PI, PI, 2 * PI});
    DomainDiscretization<Vec3d> param_domain(bs);
    domain.fill(gsf, param_domain, three_sphere_r, three_sphere_jacobian_matrix, h);

    // Calculate distances.
    KDTree<Vec<double, 4> > tree(domain.positions());
    Range<double> distances;
    double avg_dist = 0;
    for (const auto& pt : domain.positions()) {
        auto d = tree.query(pt, 10).second;
        double val = 0.0;
        for (int i = 1; i < d.size(); i++) val += sqrt(d[i]);
        val /= (d.size() - 1);

        distances.push_back(val);
        avg_dist += distances[distances.size() - 1];
    }
    avg_dist /= domain.positions().size();

    auto err_avg = std::abs(avg_dist - h) / h;
    EXPECT_LT(err_avg, 0.2);

    auto err_min = std::abs(*std::min_element(distances.begin() + 1, distances.end() - 1) - h) / h;
    EXPECT_LT(err_min, 0.2);
}

TEST(DomainEngines, GeneralSurfaceFill4DFunctionFill) {
    auto three_sphere_r = [](Vec3d t){
        double R = 10.0;
        return Vec<double, 4>(R * cos(t(0)), R * sin(t(0)) * cos(t(1)),
                              R * sin(t(0)) * sin(t(1)) * cos(t(2)),
                              R * sin(t(0)) * sin(t(1)) * sin(t(2)));
    };

    auto three_sphere_jacobian_matrix = [](Vec3d t){
        double R = 10.0;
        Eigen::Matrix<double, 4, 3> jm;
        jm.col(0) << - R * sin(t(0)), R * cos(t(0)) * cos(t(1)),
                R * cos(t(0)) * sin(t(1)) * cos(t(2)),
                R * cos(t(0)) * sin(t(1)) * sin(t(2));
        jm.col(1) << 0, - R * sin(t(0)) * sin(t(1)),
                R * sin(t(0)) * cos(t(1)) * cos(t(2)),
                R * sin(t(0)) * cos(t(1)) * sin(t(2));
        jm.col(2) << 0, 0, - R * sin(t(0)) * sin(t(1)) * sin(t(2)),
                R * sin(t(0)) * sin(t(1)) * cos(t(2));
        return jm;
    };

    auto gradient_h = [](Vec<double, 4> p){
        double h_0 = 10.0;
        return h_0 + (p[0] + p[1] + p[2] + p[3]) / 80.0;
    };

    // Define domain.
    UnknownShape<Vec<double, 4> > shape;
    DomainDiscretization<Vec<double, 4> > domain(shape);

    // Fill domain.
    GeneralSurfaceFill<Vec<double, 4>, Vec3d> gsf; gsf.seed(0);  // deterministic
    BoxShape<Vec3d> bs(Vec3d{0.0, 0.0, 0.0}, Vec3d{PI, PI, 2 * PI});
    DomainDiscretization<Vec3d> param_domain(bs);
    domain.fill(gsf, param_domain, three_sphere_r, three_sphere_jacobian_matrix, gradient_h);

    // Calculate distances.
    KDTree<Vec<double, 4> > tree(domain.positions());
    auto positions = domain.positions();
    for (int i = 1; i < positions.size() - 1 ; i++) {
        auto pt = positions[i];
        auto d = tree.query(pt, 9).second;
        double val = 0.0;
        for (int i = 1; i < d.size(); i++) val += sqrt(d[i]);
        val /= (d.size() - 1);

        double h = gradient_h(pt);
        EXPECT_GE(val, 0.75 * h);
    }
}

TEST(DomainEngines, GeneralSurfaceFillMaxPoints) {
    // Define surface.
    auto example_r = [](Vec1d t) {
        t(0) = -t(0);
        double r = pow(abs(cos(1.5 * t(0))), sin(3 * t(0)));
        return Vec2d(r * cos(t(0)), r * sin(t(0)));
    };

    auto der_example_r = [](Vec1d t) {
        t(0) = -t(0);
        double r = pow(abs(cos(1.5 * t(0))), sin(3 * t(0)));
        double der_r = (-1.5 * pow(abs(cos(1.5 * t(0))),
                                   sin(3 * t(0))) * sin(3 * t(0)) * sin(1.5 * t(0))  +
                        3 * pow(abs(cos(1.5 * t(0))),
                                sin(3 * t(0))) * cos(3 * t(0)) * cos(1.5 * t(0))
                        * log(abs(cos(1.5 * t(0))))) / cos(1.5 * t(0));

        Eigen::Matrix<double, 2, 1> jm;
        jm.col(0) << -(der_r * cos(t(0)) - r * sin(t(0))),
                -(der_r * sin(t(0)) + r * cos(t(0)));

        return jm;
    };

    double h = 0.02;

    UnknownShape<Vec2d> shape;
    DomainDiscretization<Vec2d> domain(shape);

    int max_points = 400, num_samples = 2;
    // Fill domain.
    GeneralSurfaceFill<Vec2d, Vec1d> gsf; gsf.seed(0);  // deterministic
    gsf.maxPoints(max_points).numSamples(num_samples);
    BoxShape<Vec1d> bs(Vec1d{0.0}, Vec1d{2 * PI});
    DomainDiscretization<Vec1d> param_d(bs);
    domain.fill(gsf, param_d, example_r, der_example_r, h);

    EXPECT_LE(domain.positions().size(), max_points + num_samples);
    EXPECT_GE(domain.positions().size(), max_points - num_samples);
}


TEST(DomainEngines, DISABLED_EarthModel) {
    std::ifstream fin("../testdata/topo.txt");

    double map[180][360];
    std::string is, s;
    double mx = -10000000, mn = 10000000;

    for (int i = 0; !fin.eof() && i < 180; i++) {
        fin >> is;
        std::istringstream ss(is);

        for (int j = 0; std::getline(ss, s, ',') && j < 360; j++) {
            map[i][j] = std::stoi(s);
            mx = std::max(mx, map[i][j]);
            mn = std::min(mn, map[i][j]);
            }
        }

        for (int i = 0; i < 180; i++) {
            for (int j = 0; j < 360; j++) {
                map[i][j] = (map[i][j] - mn) / (mx - mn);
        }
    }

    auto sphere_r = [](Vec2d t) {
        double R = 1.001;
        return Vec3d(R * cos(t(0)) * sin(t(1)),
                R * sin(t(0)) * sin(t(1)), R * cos(t(1)));
    };

    auto sphere_jacobian_matrix = [](Vec2d t) {
        double R = 1.001;
        Eigen::Matrix<double, 3, 2> jm;
        jm.col(0) << - R * sin(t(0)) * sin(t(1)),
                        R * cos(t(0)) * sin(t(1)), 0;
        jm.col(1) << R * cos(t(0)) * cos(t(1)),
                        R * sin(t(0)) * cos(t(1)), - R * sin(t(1));
        return jm;
    };

    auto earth_h = [&map](Vec3d p) {
        double R = p.norm();
        double theta = (90.0 - acos(p(2) / R) * 180.0 / 3.14);
        double phi = atan2(-p(0), p(1)) * 180.0 / 3.14;

        double h_0 = 0.005;
        return h_0 / (1 * map[static_cast<int>(theta - 0.5) + 90]
                             [static_cast<int>(phi - 0.5) + 180]);
    };

    UnknownShape<Vec3d> shape;
    DomainDiscretization<Vec3d> domain(shape);

    GeneralSurfaceFill<Vec3d, Vec2d> gsf; gsf.seed(0);  // deterministic
    BoxShape<Vec2d> bs({0.0, 0.0}, {2 * PI, PI});
    DomainDiscretization<Vec2d> param_domain(bs);

    domain.fill(gsf, param_domain, sphere_r, sphere_jacobian_matrix, earth_h);
}

}  // namespace mm
