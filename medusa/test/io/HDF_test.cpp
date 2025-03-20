#include <medusa/bits/io/HDF.hpp>
#include <medusa/bits/io/HDF_Eigen.hpp>
#include <Eigen/Core>
#include <medusa/bits/utils/Timer.hpp>
#include <Eigen/Sparse>
#include <medusa/bits/types/Vec.hpp>
#include <medusa/bits/types/Range.hpp>
#include <medusa/bits/domains/DomainDiscretization.hpp>
#include <medusa/bits/io/XML.hpp>

#include "gtest/gtest.h"

namespace mm {

TEST(IO, HDF5Read) {
    // open file
    std::string file_name = "test/testdata/test_hdf5_read.h5";
    HDF reader(file_name, HDF::READONLY);
    EXPECT_TRUE(reader.isFileOpen());
    EXPECT_TRUE(reader.isGroupOpen());

    // read parameters in /what
    reader.openGroup("/group1");
    EXPECT_EQ("/group1", reader.groupName());
    reader.openGroup("how");
    EXPECT_EQ("/group1/how", reader.groupName());
    reader.openGroup("/");
    EXPECT_EQ("/", reader.groupName());
    reader.openGroup("group1");
    EXPECT_EQ("/group1", reader.groupName());

    reader.close();
    reader.open(file_name, "/what", HDF::READONLY);

    EXPECT_EQ("/what", reader.groupName());
    EXPECT_EQ(file_name, reader.filename());

    std::string read_str = reader.readStringAttribute("object");
    EXPECT_EQ(std::string("PVOL\0", 5), read_str);
    std::string date_str = reader.readNullTerminatedStringAttribute("date");
    EXPECT_EQ("20150501", date_str);
    std::string time_str = reader.readNullTerminatedStringAttribute("time");
    EXPECT_EQ("000000", time_str);
    std::string source_str = reader.readNullTerminatedStringAttribute("source");
    EXPECT_EQ("WMO:14024,RAD:SI41,PLC:Lisca,NOD:silis", source_str);

    reader.openGroup("/group1/data1");
    std::vector<std::string> groups = {"how", "what"}, datasets = {"data"};
    EXPECT_EQ(groups, reader.groups());
    EXPECT_EQ(datasets, reader.datasets());

    reader.openGroup("/where");
    double height = reader.readDoubleAttribute("height");
    EXPECT_DOUBLE_EQ(950.0, height);
    double lat = reader.readDoubleAttribute("lat");
    EXPECT_DOUBLE_EQ(46.06776997447014, lat);
    double lon = reader.readDoubleAttribute("lon");
    EXPECT_DOUBLE_EQ(15.28489999473095, lon);

    reader.openGroup("/group1/how");
    double avgpwr = reader.readDoubleAttribute("avgpwr");
    EXPECT_DOUBLE_EQ(135.0, avgpwr);

    reader.openGroup("/group2/how");
    double avgpwr2 = reader.readDoubleAttribute("avgpwr");
    EXPECT_DOUBLE_EQ(135.0, avgpwr2);

    reader.openGroup("/group1/data1/what");
    std::string quantity = reader.readNullTerminatedStringAttribute("quantity");
    EXPECT_EQ("DBZH", quantity);
    reader.openGroup("/group1/data2/what");
    quantity = reader.readNullTerminatedStringAttribute("quantity");
    EXPECT_EQ("VRAD", quantity);
}

class HDF5WriteTest : public ::testing::Test {
  public:
    HDF file;
    static const char file_name[];
  protected:
    void SetUp() override {
        file.openFile(file_name, HDF::DESTROY);
        file.openGroup("/test");
    }

    void TearDown() override {
        file.closeGroup();
        file.closeFile();
        std::remove(file_name);
    }
};
const char HDF5WriteTest::file_name[] = "test/testdata/test_hdf5_write.h5";

TEST_F(HDF5WriteTest, Reopen) {
    file.openGroup("/test/test2/test3");
    double attr = 1.2;
    file.writeDoubleAttribute("a", attr);
    file.close();
    EXPECT_FALSE(file.isFileOpen());
    EXPECT_FALSE(file.isGroupOpen());
    file.reopenFile();
    EXPECT_TRUE(file.isFileOpen());
    EXPECT_FALSE(file.isGroupOpen());
    file.reopenGroup();
    EXPECT_TRUE(file.isFileOpen());
    EXPECT_TRUE(file.isGroupOpen());
    EXPECT_EQ(attr, file.readDoubleAttribute("a"));
}

TEST_F(HDF5WriteTest, Attribute) {
    // root folder
    file.openGroup("/");
    file.writeIntAttribute("int_attribute", 5);
    EXPECT_EQ(5, file.readIntAttribute("int_attribute"));

    // ordinary folder
    file.openGroup("/test");
    // INT
    file.writeIntAttribute("int_attribute", -10);
    EXPECT_EQ(-10, file.readIntAttribute("int_attribute"));

    // folder within existing folder
    file.openGroup("/test/test2");
    // DOUBLE
    file.writeDoubleAttribute("double_attribute", 5.1);
    EXPECT_DOUBLE_EQ(5.1, file.readDoubleAttribute("double_attribute"));
    file.writeDoubleAttribute("double_attribute", 10.1245, true);
    EXPECT_DOUBLE_EQ(10.1245, file.readDoubleAttribute("double_attribute"));

    // new double folder
    file.openGroup("/test1/test1");
    // FLOAT
    file.writeFloatAttribute("float_attribute", 5.1f);
    EXPECT_FLOAT_EQ(5.1f, file.readFloatAttribute("float_attribute"));
    file.writeFloatAttribute("float_attribute", 10.1245f, true);
    EXPECT_FLOAT_EQ(10.1245f, file.readFloatAttribute("float_attribute"));

    // new double folder
    file.openGroup("/test1/test1");
    // BOOL
    file.writeBoolAttribute("bool_attribute", true);
    EXPECT_TRUE(file.readBoolAttribute("bool_attribute"));
    file.writeBoolAttribute("bool_attribute", false, true);
    EXPECT_FALSE(file.readBoolAttribute("bool_attribute"));

    // existing folder
    file.openGroup("/test/test2");
    // STRING
    std::string s = "I have a null terminator in the middle";
    s[8] = '\0';
    file.writeStringAttribute("tricky_string", s);
    EXPECT_EQ(s, file.readStringAttribute("tricky_string"));
    file.writeStringAttribute("string_attribute", "Test string");
    EXPECT_EQ("Test string", file.readStringAttribute("string_attribute"));
    file.writeStringAttribute("string_attribute", "Second test string 1234.", true);
    EXPECT_EQ("Second test string 1234.", file.readStringAttribute("string_attribute"));
}

TEST_F(HDF5WriteTest, DoubleArray) {
    std::vector<double> ret, data = {1.0, 10.234, 2.4};
    file.writeDoubleArray("vector_doubles", data);
    ret = file.readDoubleArray("vector_doubles");
    ASSERT_EQ(data.size(), ret.size());
    for (size_t i = 0; i < ret.size(); i++) {
        EXPECT_DOUBLE_EQ(data[i], ret[i]);
    }

    data = {5.0, 1234.234, 556};
    file.writeDoubleArray("vector_doubles", data, true);
    ret = file.readDoubleArray("vector_doubles");
    ASSERT_EQ(data.size(), ret.size());
    for (size_t i = 0; i < ret.size(); i++) {
        EXPECT_DOUBLE_EQ(data[i], ret[i]);
    }
}

TEST_F(HDF5WriteTest, DoubleArray2D) {
    std::vector<std::vector<double>> data2d = {
            {0.12, 3.14157, 12.7}, {2.12, 2.6, 18.65}, {3.12, 4.6, 99.65}, {5, 5, 5}};
    file.writeDouble2DArray("vec_vec_double", data2d);
    std::vector<std::vector<double>> ret2d = file.readDouble2DArray("vec_vec_double");
    ASSERT_EQ(data2d.size(), ret2d.size());
    for (size_t i = 0; i < data2d.size(); i++) {
        ASSERT_EQ(data2d[i].size(), ret2d[i].size());
        for (size_t j = 0; j < data2d[i].size(); j++) {
            EXPECT_EQ(data2d[i][j], ret2d[i][j]);
        }
    }
}

TEST_F(HDF5WriteTest, DoubleArray3D) {
    std::vector<std::vector<std::vector<double>>> data3d = {
            {{0.12, 3.14157, 12.7}, {2.12, 2.6, 18.65}, {3.12, 4.6, 99.65}, {5, 5, 5}},
            {{-0.12, -3.14157, -12.7}, {-2.12, -2.6, -18.65}, {-3.12, -4.6, -99.65}, {-5, -5, -5}}};
    file.writeDouble3DArray("v3_double", data3d);
    std::vector<std::vector<std::vector<double>>> ret3d = file.readDouble3DArray("v3_double");
    ASSERT_EQ(data3d.size(), ret3d.size());
    for (size_t i = 0; i < data3d.size(); i++) {
        ASSERT_EQ(data3d[i].size(), ret3d[i].size());
        for (size_t j = 0; j < data3d[i].size(); j++) {
            ASSERT_EQ(data3d[i][j].size(), ret3d[i][j].size());
            for (size_t k = 0; k < data3d[i][j].size(); k++) {
                EXPECT_EQ(data3d[i][j][k], ret3d[i][j][k]);
            }
        }
    }
}

TEST_F(HDF5WriteTest, FloatArray) {
    std::vector<float> ret;
    std::vector<float> data = {1.0, 10.234, 2.4};
    file.writeFloatArray("vector_floats", data);
    ret = file.readFloatArray("vector_floats");
    ASSERT_EQ(data.size(), ret.size());
    for (size_t i = 0; i < ret.size(); i++) {
        EXPECT_EQ(data[i], ret[i]);
    }

    data = {5.0, 1234.234, 556};
    file.writeFloatArray("vector_floats", data, true);
    ret = file.readFloatArray("vector_floats");
    ASSERT_EQ(data.size(), ret.size());
    for (size_t i = 0; i < ret.size(); i++) {
        EXPECT_EQ(data[i], ret[i]);
    }
}

TEST_F(HDF5WriteTest, FloatArray2D) {
    std::vector<std::vector<float>> data2d = {
            {0.12, 3.14157, 12.7}, {2.12, 2.6, 18.65}, {3.12, 4.6, 99.65}, {5, 5, 5}};
    file.writeFloat2DArray("vec_vec_float", data2d);
    std::vector<std::vector<float>> ret2d = file.readFloat2DArray("vec_vec_float");
    ASSERT_EQ(data2d.size(), ret2d.size());
    for (size_t i = 0; i < data2d.size(); i++) {
        ASSERT_EQ(data2d[i].size(), ret2d[i].size());
        for (size_t j = 0; j < data2d[i].size(); j++) {
            EXPECT_EQ(data2d[i][j], ret2d[i][j]);
        }
    }
}

TEST_F(HDF5WriteTest, FloatArray3D) {
    std::vector<std::vector<std::vector<float>>> data3d = {
            {{0.12, 3.14157, 12.7}, {2.12, 2.6, 18.65}, {3.12, 4.6, 99.65}, {5, 5, 5}},
            {{-0.12, -3.14157, -12.7}, {-2.12, -2.6, -18.65}, {-3.12, -4.6, -99.65}, {-5, -5, -5}}};
    file.writeFloat3DArray("v3_float", data3d);
    std::vector<std::vector<std::vector<float>>> ret3d = file.readFloat3DArray("v3_float");
    ASSERT_EQ(data3d.size(), ret3d.size());
    for (size_t i = 0; i < data3d.size(); i++) {
        ASSERT_EQ(data3d[i].size(), ret3d[i].size());
        for (size_t j = 0; j < data3d[i].size(); j++) {
            ASSERT_EQ(data3d[i][j].size(), ret3d[i][j].size());
            for (size_t k = 0; k < data3d[i][j].size(); k++) {
                EXPECT_EQ(data3d[i][j][k], ret3d[i][j][k]);
            }
        }
    }
}

TEST_F(HDF5WriteTest, IntArray) {
    std::vector<int> ret;
    std::vector<int> data = {1, 10, 2};
    file.writeIntArray("vector_ints", data);
    ret = file.readIntArray("vector_ints");
    ASSERT_EQ(data.size(), ret.size());
    for (size_t i = 0; i < ret.size(); i++) {
        EXPECT_EQ(data[i], ret[i]);
    }
}

TEST_F(HDF5WriteTest, IntArray2D) {
    std::vector<std::vector<int>> data2d = {{12, 3, 12}, {2, 2, 18}, {3, 4, 99}, {5, 5, 5}};
    file.writeInt2DArray("vec_vec_int", data2d);
    std::vector<std::vector<int>> ret2d = file.readInt2DArray("vec_vec_int");
    ASSERT_EQ(data2d.size(), ret2d.size());
    for (size_t i = 0; i < data2d.size(); i++) {
        ASSERT_EQ(data2d[i].size(), ret2d[i].size());
        for (size_t j = 0; j < data2d[i].size(); j++) {
            EXPECT_EQ(data2d[i][j], ret2d[i][j]);
        }
    }
}

TEST_F(HDF5WriteTest, IntArray3D) {
    std::vector<std::vector<std::vector<int>>> data3d = {
            {{12, 3, 12}, {2, 2, 18}, {3, 4, 99}, {5, 5, 5}},
            {{12, 3, 12}, {2, 2, 18}, {3, 4, 99}, {5, 5, 5}}};
    file.writeInt3DArray("v3_int", data3d);
    std::vector<std::vector<std::vector<int>>> ret3d = file.readInt3DArray("v3_int");
    ASSERT_EQ(data3d.size(), ret3d.size());
    for (size_t i = 0; i < data3d.size(); i++) {
        ASSERT_EQ(data3d[i].size(), ret3d[i].size());
        for (size_t j = 0; j < data3d[i].size(); j++) {
            ASSERT_EQ(data3d[i][j].size(), ret3d[i][j].size());
            for (size_t k = 0; k < data3d[i][j].size(); k++) {
                EXPECT_EQ(data3d[i][j][k], ret3d[i][j][k]);
            }
        }
    }
}

TEST_F(HDF5WriteTest, ArrayWithVec) {
    Range<Vec<double, 3>> vector_data;
    for (auto i = 1; i < 100; ++i) vector_data.push_back(Vec<double, 3>{0.1, 0.2, 0.3});
    file.writeDouble2DArray("vector_data", vector_data);

    std::vector<std::vector<double>> tt = file.readDouble2DArray("vector_data");
    for (auto i = 0; i < vector_data.size(); ++i) {
        for (auto j = 0; j < vector_data[0].size(); ++j)
            EXPECT_DOUBLE_EQ(tt[i][j], vector_data[i][j]);
    }
}

TEST_F(HDF5WriteTest, Array3D) {
    Range<Range<Vec<double, 3>>> data(5, Range<Vec<double, 3>>(4));
    for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 4; ++j) {
    for (int k = 0; k < 3; ++k) {
        data[i][j][k] = 0.5*i - 2.1*j + 1.2*k;
    }}}
    file.write3DArray<double>("data3d", data, H5T_NATIVE_DOUBLE, false);
}

TEST_F(HDF5WriteTest, Eigen) {
    Eigen::Matrix2d M; M.setRandom();

    file.writeEigen("M", M);
    Eigen::MatrixXd M2 = file.readEigen("M");
    EXPECT_TRUE((M-M2).isZero(0));

    file.writeEigen("M", M*M, true);
    M2 = file.readEigen("M");
    EXPECT_TRUE((M*M - M2).isZero(0));

    Eigen::MatrixX2d v(12, 2);
    v.setRandom();
    file.writeEigen("testEigen", v, true);
    auto v2 = file.readDouble2DArray("testEigen");

    for (int i = 0; i < v.rows(); ++i) {
        ASSERT_EQ(static_cast<int>(v2[i].size()), v.cols());
        for (int j = 0; j < v.cols(); ++j) {
            EXPECT_EQ(v(i, j), v2[i][j]);
        }
    }

    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> M3(3, 2);
    M3.setRandom();
    file.writeEigen("M3", M3);
    auto M4 = file.readEigen<float>("M3");
    EXPECT_TRUE((M3 - M4).isZero(0));

    int m = 5, n = 3;
    std::vector<std::vector<float>> M5(m);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            M5[i].push_back(2.32f*j+0.23f*i);
        }
    }
    file.writeFloat2DArray("M5", M5);
    auto M6 = file.readEigen<float>("M5");

    ASSERT_EQ(m, M6.rows());
    ASSERT_EQ(n, M6.cols());
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            EXPECT_EQ(M6(i, j), M5[i][j]);
        }
    }
}

TEST_F(HDF5WriteTest, SparseMatrix) {
    Eigen::SparseMatrix<double> M(10, 10);
    M.coeffRef(0, 0) = 4.3;
    M.coeffRef(4, 5) = -2.1;
    file.writeSparseMatrix("M1", M);
    file.writeSparseMatrix("M0", M, false);

    auto M1 = file.readDouble2DArray("M1");
    std::vector<std::vector<double>> M1_expected = {{1, 1, 4.3}, {5, 6, -2.1}};
    EXPECT_EQ(M1_expected, M1);
    auto M0 = file.readDouble2DArray("M0");
    std::vector<std::vector<double>> M0_expected = {{0, 0, 4.3}, {4, 5, -2.1}};
    EXPECT_EQ(M0_expected, M0);
}

TEST_F(HDF5WriteTest, Timer) {
    Timer t;
    t.addCheckPoint("a");
    t.addCheckPoint("b");
    t.addCheckPoint("c");
    file.writeTimer("time", t);
    file.openGroup("/test/time");
    double total = file.readDoubleAttribute("total");
    EXPECT_LE(total, 1e-2);
    EXPECT_GE(total, 0);

    double ab = file.readDoubleAttribute("a-b");
    EXPECT_LE(ab, 1e-2);
    EXPECT_GE(ab, 0);

    double bc = file.readDoubleAttribute("b-c");
    EXPECT_LE(bc, 1e-2);
    EXPECT_GE(bc, 0);
}

TEST_F(HDF5WriteTest, Conf) {
    XML conf("test/testdata/test_conf.xml");
    file.writeXML("conf", conf);

    file.openGroup("/test/conf");
    double actual = file.readDoubleAttribute("mls.m");
    EXPECT_EQ(3, actual);
    actual = file.readDoubleAttribute("solver.droptol");
    EXPECT_EQ(1e-5, actual);
    std::string data = file.readStringAttribute("mls.basis_type");
    EXPECT_EQ("mon", data);
    data = file.readStringAttribute("num.nxs");
    EXPECT_EQ("10, 12, 14, 16", data);
    data = file.readStringAttribute("meta.filename");
    EXPECT_EQ("one_domain_implicit_smooth.h5", data);
}

TEST_F(HDF5WriteTest, Domain) {
    UnknownShape<Vec2d> shape;
    DomainDiscretization<Vec2d> d(shape);
    d.addInternalNode({1.2, -0.3}, 1);
    d.addInternalNode({0.4, 0.0}, 2);
    d.addBoundaryNode({0.5, -0.1}, -1, {-1, 0});
    d.addInternalNode({12.3, 1.8}, 3);
    d.addBoundaryNode({4.5, 3.1}, -2, {0, 1});

    file.openGroup("/test/domain_test");
    file.writeDomain("domain", d);

    file.openGroup("domain");
    EXPECT_EQ(5, file.readIntAttribute("N"));
    std::vector<int> expected_types = {1, 2, -1, 3, -2};
    EXPECT_EQ(expected_types, file.readIntArray("types"));
    Eigen::MatrixXd M(5, 2); M << 1.2, -0.3, 0.4, 0.0, 0.5, -0.1, 12.3, 1.8, 4.5, 3.1;
    EXPECT_TRUE((M - file.readEigen("pos")).isZero(0));
    Eigen::MatrixXd normals(2, 2); normals << -1, 0, 0, 1;
    EXPECT_TRUE((normals - file.readEigen("normals")).isZero(0));
    std::vector<int> bmap = {-1, -1, 0, -1, 1};
    EXPECT_EQ(bmap, file.readIntArray("bmap"));
}

TEST_F(HDF5WriteTest, Atomic) {
    file.close();

    file.setGroupName("asdf");
    EXPECT_EQ("/asdf", file.groupName());
    file.setGroupName("asdf/");
    EXPECT_EQ("/asdf", file.groupName());
    file.setGroupName("/");
    EXPECT_EQ("/", file.groupName());

    EXPECT_FALSE(file.isGroupOpen());
    EXPECT_FALSE(file.isFileOpen());
    file.atomic().writeDoubleAttribute("a", 12.3);
    EXPECT_FALSE(file.isGroupOpen());
    EXPECT_FALSE(file.isFileOpen());
    double a = file.atomic().readDoubleAttribute("a");
    EXPECT_FALSE(file.isGroupOpen());
    EXPECT_FALSE(file.isFileOpen());
    EXPECT_EQ(a, 12.3);
}

TEST(IO, HDFDeath) {
    HDF reader;
    EXPECT_DEATH(reader.openFile("sadfasdf", HDF::READONLY), "must exist and be accessible");
    reader.openFile("test/testdata/test_hdf5_read.h5");
    reader.openGroup("/what");
    EXPECT_DEATH(reader.writeStringAttribute("object", "content"), "already exists");
}

TEST(IO, HDFUsageExample) {
    if (std::ifstream("test/testdata/new_file.h5").good()) {
        std::remove("test/testdata/new_file.h5");
    }
    /// [HDF usage example]
    // open new empty file and group '/' automatically
    HDF hdf("test/testdata/filename.h5", HDF::DESTROY);
    hdf.writeStringAttribute("text", "This is the content");  // written to group '/'

    hdf.openGroup("group1");
    std::vector<double> a = {1.2, -4.5, 4.5};
    hdf.writeDoubleArray("a", a);
    hdf.close();

    hdf.reopen();  // reopen previously closed file and group
    hdf.openGroup("/");
    std::string s = hdf.readStringAttribute("text");
    hdf.writeStringAttribute("text", "New content", true);  // enable overwrite

    hdf.setFilename("test/testdata/new_file.h5");
    hdf.setGroupName("group_name");  // set current group name and filename, but do not open yet

    std::vector<std::vector<double>> b(5, a);  // written as dataset of size 5x3
    hdf.atomic().writeFloat2DArray("b", b);  // opens, writes and closes, the value is cast to float

    Eigen::MatrixXd M(5, 3);
    M << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15;
    hdf.atomic().writeEigen("M", M);
    auto M2 = hdf.atomic().readEigen("M");  // opens, reads, closes

    std::cout << hdf << std::endl;
    // Timer, DomainDiscretization and XML classes can also be written directly
    /// [HDF usage example]
    s = "";
    std::remove("test/testdata/new_file.h5");
    EXPECT_TRUE((M2 - M).isZero(0));
}

}  // namespace mm
