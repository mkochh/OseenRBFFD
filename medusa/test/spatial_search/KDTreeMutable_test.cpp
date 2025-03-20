#include <medusa/bits/spatial_search/KDTreeMutable_fwd.hpp>

#include "gtest/gtest.h"

#include <medusa/bits/types/Vec.hpp>

namespace mm {

TEST(KDTree, DISABLED_MutableUsageExample) {
    ///  [KDTreeMutable usage example]
    Range<Vec3d> points = {{5, 4, 0}, {4, 2, 1}, {7, 6, 9}, {2, 2, 1},
                           {8, 0, 5}, {5, 7, 0}, {3, 3, 8}, {9, 7, 3}, {2, 2, 6}, {2, 0, 6}};
    KDTreeMutable<Vec3d> tree(points);  //  Build kd-tree
    Range<double> distances;
    Range<int> closest;
    //  Get 3 points closest to {4.1, 3, 0.5}
    std::tie(closest, distances) = tree.query({4.1, 3, 0.5}, 3);
    //  Insertion and removal
    Vec3d new_point = {2, 30, 45};
    tree.insert(new_point);
    tree.remove(2);  // remove by index (assigned at insertion)
    std::cout << tree << std::endl;
    ///  [KDTreeMutable usage example]
}

TEST(KDTree, Mutable3DQuery) {
    Range<Vec3d> points = {
            {5, 4, 0}, {4, 2, 1}, {7, 6, 9}, {2, 2, 1}, {8, 0, 5}, {5, 7, 0},
            {3, 3, 8}, {9, 7, 3}, {2, 2, 6}, {2, 0, 6}};
    KDTreeMutable<Vec3d> kd3(points);
    EXPECT_EQ(points.size(), kd3.size());
    EXPECT_EQ(std::vector<int>({1, 0}), kd3.query({4.1, 3, 0.5}, 2).first);
    EXPECT_EQ(std::vector<int>({3, 1, 8}), kd3.query({2, 2, 3}, 3).first);
}


TEST(KDTree, Reset) {
    Range<Vec3d> points = {
            {5, 4, 0}, {4, 2, 1}, {7, 6, 9}, {2, 2, 1}, {8, 0, 5}, {5, 7, 0},
            {3, 3, 8}, {9, 7, 3}, {2, 2, 6}, {2, 0, 6}};
    KDTreeMutable<Vec3d> kd3(points);
    EXPECT_EQ(points.size(), kd3.size());
    EXPECT_EQ(std::vector<int>({1, 0}), kd3.query({4.1, 3, 0.5}, 2).first);
    EXPECT_EQ(std::vector<int>({3, 1, 8}), kd3.query({2, 2, 3}, 3).first);
    kd3.reset({{1, 2, 3}, {4, 5, 6}});
    EXPECT_EQ(std::vector<int>({0, 1}), kd3.query({-1, -2, -3}, 2).first);
    EXPECT_EQ(std::vector<double>({1*1+2*2+3*3, 4*4+5*5+6*6}), kd3.query({0, 0, 0}, 2).second);
}

TEST(KDTree, Mutable2DQuery) {
    // Test 2D with 50 points
    Range<Vec2d> points = {
            {-0.6441209666616923, 0.5850353695296744}, {-0.6201378980490286, 0.7811830957022021},
            {-0.6784344365993658, 0.9839667035596042}, {0.044496636767509035, 0.4377470158269272},
            {-0.5274564029904492, 0.5696927123044102}, {-0.6259567688687637, -0.7818614711722951},
            {0.9162830614939184, 0.6489213976078669}, {0.7591192942947527, -0.3433314817836943},
            {-0.23782320760874454, 0.1898945196158186}, {-0.5441698103664361, 0.006385132320767317},
            {-0.3654415676489693, 0.06746853328508307}, {0.931189791692042, 0.9181001537912621},
            {-0.3427312891982608, -0.6826010311723985}, {0.891723271821818, -0.8508989951879877},
            {0.11394616056891804, 0.84630479540056}, {-0.5223541537898675, -0.026438581866889965},
            {0.425710227280224, 0.232320814224092}, {-0.12153165386201059, -0.09554948988436873},
            {-0.5360668060424465, 0.8346421613694428}, {-0.48815760323401425, -0.47632741760208863},
            {-0.7355755648085145, 0.021559400675411178}, {0.16684622124842763, -0.5180297967184262},
            {0.2750341401977108, 0.14027191172883025}, {0.6550291085853586, -0.0906298270828052},
            {0.11672354315833089, 0.6864647718719339}, {0.4319012792644372, -0.5056121884571734},
            {-0.3111112525746833, 0.4841999490074793}, {0.9537933710310402, -0.4121810950684375},
            {-0.13430435085483827, -0.7063126623794711}, {0.8509710074195942, 0.4342793339962312},
            {0.8387012530238107, -0.8797447958759539}, {-0.788610725923649, -0.28402465290256296},
            {0.3781644510827291, -0.19185822357265714}, {0.3675280148118516, -0.19959444568922202},
            {0.5441631640348303, 0.9923042351466602}, {0.4379307570126687, -0.7976305308639227},
            {0.6192156716466921, 0.5095818625700927}, {0.4538538344127916, 0.24132983301916466},
            {-0.920957301103486, 0.7946521172232996}, {-0.4620339775822764, -0.6725241136453264},
            {0.6164340529959673, -0.370457261186975}, {0.6353232887100764, -0.4527336246027056},
            {0.5129962829311199, 0.2576340035951956}, {-0.3355508932427227, 0.20820573070479842},
            {-0.03557537269552791, 0.801250103966441}, {-0.6249260602423057, 0.01901038540320088},
            {0.9052024446464844, -0.37370939241714884}, {-0.08417808324593734, -0.4918951358215899},
            {-0.7408521674604813, -0.38351394950328266}, {0.8804020099004968, 0.2246660525735673}};
    int rs;
    KDTreeMutable<Vec2d> kd2(points);

    EXPECT_EQ(points.size(), kd2.size());
    //  EXPECT_EQ(2u, kd2.query(points[33], 0.05, 3).first.size());
    //  EXPECT_EQ(2u, kd2.query(points[33], 0.05, 2).first.size());
    auto result = kd2.query(points[33], 3);
    EXPECT_EQ(std::vector<int>({33, 32, 40}), result.first);
    rs = result.first.size();
    for (int i = 0; i < rs; i++) {
        EXPECT_DOUBLE_EQ(
                (points[33] - points[result.first[i]]).squaredNorm(), result.second[i]);
        EXPECT_LE(result.second[i], result.second[rs-1]);
    }

    //  EXPECT_EQ(3u, kd2.query(points[17], 0.09, 3).first.size());
    //  EXPECT_EQ(2u, kd2.query(points[17], 0.09, 2).first.size());
    result = kd2.query(points[17], 3);
    EXPECT_EQ(std::vector<int>({17, 10, 8}), result.first);
    rs = result.first.size();
    for (int i = 0; i < rs; i++) {
        EXPECT_DOUBLE_EQ(
                (points[17] - points[result.first[i]]).squaredNorm(), result.second[i]);
        EXPECT_LE(result.second[i], result.second[rs-1]);
    }
    //  EXPECT_EQ(4u, kd2.query({0., 0}, 0.14, 5).first.size());
    //  EXPECT_EQ(2u, kd2.query({0., 0}, 0.14, 2).first.size());
    result = kd2.query({0., 0}, 5);
    EXPECT_EQ(std::vector<int>({17, 8, 22, 10, 43}), result.first);
    rs = result.first.size();
    for (int i = 0; i < rs; i++) {
        EXPECT_DOUBLE_EQ(
                (Vec2d({0, 0}) - points[result.first[i]]).squaredNorm(), result.second[i]);
        EXPECT_LE(result.second[i], result.second[rs-1]);
    }

    //  EXPECT_EQ(3u, kd2.query({0.7, 0.7}, 0.101, 5).first.size());
    //  EXPECT_EQ(2u, kd2.query({0.7, 0.7}, 0.101, 2).first.size());
    result = kd2.query({0.7, 0.7}, 5);
    EXPECT_EQ(std::vector<int>({36, 6, 29, 11, 34}), result.first);
    rs = result.first.size();
    for (int i = 0; i < rs; i++) {
        EXPECT_DOUBLE_EQ(
                (Vec2d({0.7, 0.7}) - points[result.first[i]]).squaredNorm(), result.second[i]);
        EXPECT_LE(result.second[i], result.second[rs-1]);
    }
}

TEST(KDTree, DISABLED_MutableDeathTests) {
    KDTreeMutable<Vec2d> kd2;
    EXPECT_DEATH(kd2.query({0., 2.0}, 3), "Not enough points in the tree, you requested 3 points, "
                                          "but the tree contains only 0 points.");
    EXPECT_DEATH(kd2.insert(Vec2d{sqrt(-1), 0}), "Invalid point.");
}

TEST(KDTree, Mutableinsert_remove) {
    Range<Vec2d> points = {{2, 3}, {1, 5}, {7, 8}, {0, -2}, {9, 12}, {4, 15}, {15, 4}};
    KDTreeMutable<Vec2d> tree(points);
    EXPECT_EQ(7, tree.size());

    Vec2d new_point({2, 30});
    tree.insert(new_point);
    EXPECT_EQ(8, tree.size());

    tree.remove(0);
    EXPECT_EQ(7, tree.size());

    tree.remove(0);  // Removing a point again is valid, but has no effect.
    EXPECT_EQ(7, tree.size());
}

TEST(KDTree, Mutable2DQuery_2) {
    Range<Vec2d> points = {{2, 4}, {2, 3}, {1, 5}, {7, 8}, {0, -2}, {9, 12}, {4, 15}, {15, 4}};
    KDTreeMutable<Vec2d> tree(points);
    Vec2d reference_point({2, 3});
    mm::Range<int> closest;
    mm::Range<double> distances;
    std::tie(closest, distances) = tree.query(reference_point, 1);
    EXPECT_EQ(1 , closest[0]);
    EXPECT_EQ(0.0, distances[0]);
}

TEST(KDTree, ExistsPointInSphere) {
    Range<Vec2d> points = {{2, 4}, {2, 3}, {1, 5}, {7, 8}, {0, -2}, {9, 12}, {4, 15}, {15, 4}};
    KDTreeMutable<Vec2d> tree;
    EXPECT_FALSE(tree.existsPointInSphere({0.0, 0.0}, 1000));
    tree.insert(points);
    Vec2d reference_point({2, 3.1});
    EXPECT_TRUE(tree.existsPointInSphere(reference_point, 0.15));
    EXPECT_FALSE(tree.existsPointInSphere(reference_point, 0.15*0.15));
}

TEST(KDTree, Mutable3DQuery_2) {
    Range<Vec3d> points = {{2, 3, 7}, {2, 4, -9.5}, {1, -5, 7.2},
                           {7, 8, -40.6}, {0, -2, 20.9}, {9, 12, -102.1}, {4, 15, 99.1},
                           {15, 4, -0.4}, { -4.2, 6, -9.8}, { -1.5, -4.0, -7.0}};
    KDTreeMutable<Vec3d> tree(points);
    Vec3d reference_point({1, 3.0, 8.0});
    mm::Range<int> closest_points;
    mm::Range<double> distances;
    std::tie(closest_points, distances) = tree.query(reference_point, 2);

    mm::Range<int>solution({0, 2});
    mm::Range<double> dist_solution({2.0,  0.8 * 0.8 + 8.0 * 8.0});
    EXPECT_EQ(solution, closest_points);
    EXPECT_EQ(dist_solution, distances);
}

}  // namespace mm
