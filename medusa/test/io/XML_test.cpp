#include <medusa/bits/io/XML.hpp>

#include "gtest/gtest.h"

namespace mm {

TEST(IO, XMLValidPath) {
    XML xml("test/testdata/test_conf.xml");
    EXPECT_DEATH(xml.exists(""), "Path name should not be empty");
    EXPECT_DEATH(xml.exists("test."), "Path .* contains an empty attribute name");
    EXPECT_DEATH(xml.exists(".test"), "Path .* contains an empty element name");
    EXPECT_DEATH(xml.exists("test.test..test"), "Path .* contains an empty element name");
}

TEST(IO, XMLExists) {
    XML xml("test/testdata/test_conf.xml");
    EXPECT_TRUE(xml.exists("h"));
    EXPECT_TRUE(xml.exists("num.nxs"));
    EXPECT_TRUE(xml.exists("solver.droptol"));
    EXPECT_TRUE(xml.exists("mls.m"));
    EXPECT_TRUE(xml.exists("mls.abc.def.m"));  // abc as an element
    EXPECT_TRUE(xml.exists("mls.abc"));  // abc as an attribute

    EXPECT_FALSE(xml.exists("solver.droptols"));
    EXPECT_FALSE(xml.exists("solve.droptols"));
    EXPECT_FALSE(xml.exists("mls.abc.a"));
    EXPECT_FALSE(xml.exists("qwer.asdf.qwer"));
    EXPECT_FALSE(xml.exists("h.h"));
    EXPECT_FALSE(xml.exists("mls.m.a"));
}

TEST(IO, XMLCopy) {
    XML xml("test/testdata/test_conf.xml");
    XML xml_copy = xml;
    EXPECT_EQ(xml.get<int>("h"), xml_copy.get<int>("h"));
    EXPECT_EQ(xml.get<std::string>("num.nxs"), xml_copy.get<std::string>("num.nxs"));
    EXPECT_EQ(xml.get<double>("solver.droptol"), xml_copy.get<double>("solver.droptol"));
    EXPECT_EQ(xml.get<int>("mls.m"), xml_copy.get<int>("mls.m"));
}

TEST(IO, XMLGet) {
    XML xml("test/testdata/test_conf.xml");
    EXPECT_EQ(6, xml.get<int>("h"));
    EXPECT_EQ("10, 12, 14, 16", xml.get<std::string>("num.nxs"));
    EXPECT_EQ(1e-5, xml.get<double>("solver.droptol"));
    EXPECT_EQ("3", xml.get<std::string>("mls.m"));
    EXPECT_EQ(3, xml.get<int>("mls.m"));
    EXPECT_EQ(3.0f, xml.get<float>("mls.m"));
    EXPECT_EQ(3.0, xml.get<double>("mls.m"));
    EXPECT_EQ(5, xml.get<int>("mls.abc.def.m"));
    EXPECT_EQ("9", xml.get<std::string>("mls.abc"));
}

TEST(IO, XMLSet) {
    XML xml("test/testdata/test_conf.xml");

    double r = 12.1372864128934;
    xml.set("r", r);
    EXPECT_EQ(r, xml.get<double>("r"));

    std::string s = "   \n   ";
    xml.set("mls.abc.s", s);
    EXPECT_EQ(s, xml.get<std::string>("mls.abc.s"));

    int i = -234;
    EXPECT_DEATH(xml.set("mls.abc.def.m", i), "Attribute on path '.*' already exists");
    xml.set("mls.abc.def.m", i, true);
    EXPECT_EQ(i, xml.get<int>("mls.abc.def.m"));

    bool b = false;
    xml.set("random.elements.should.be.created", b);
    EXPECT_EQ(b, xml.get<bool>("random.elements.should.be.created"));
    xml.set("random.elements.should.be.created", !b, true);
    EXPECT_EQ(!b, xml.get<bool>("random.elements.should.be.created"));
}

TEST(IO, XMLGetAll) {
    XML xml("test/testdata/test_conf.xml");

    double r = 12.1372864128934;
    xml.set("r", r);
    EXPECT_EQ(r, xml.get<double>("r"));
    std::vector<std::pair<std::string, std::string>> expected = {
            {"h",              "6"},
            {"r",              "12.1372864128934"},
            {"mls.m",          "3"},
            {"mls.basis_type", "mon"},
            {"mls.abc",        "9"},
            {"mls.abc.def.a",  "5"},
            {"mls.abc.def.m",  "5"},
            {"mls.abc.def2.a", "6"},
            {"solver.droptol", "1e-5"},
            {"num.nxs",        "10, 12, 14, 16"},
            {"meta.filename",  "one_domain_implicit_smooth.h5"}};
    auto actual = xml.getAll();
    std::sort(expected.begin(), expected.end());
    std::sort(actual.begin(), actual.end());
    EXPECT_EQ(expected, actual);
}

TEST(IO, XMLDoubleLoad) {
    XML xml("test/testdata/test_conf.xml");
    xml.set("r", 23.4);
    EXPECT_EQ(23.4, xml.get<double>("r"));

    xml.load("test/testdata/test_conf.xml");
    EXPECT_FALSE(xml.exists("r"));
}

TEST(IO, DISABLED_XMLUsageExample) {
    /// [XML usage example]
    XML xml("test/testdata/test_conf.xml");
    if (xml.exists("element1.attr")) {
        int a = xml.get<int>("element1.attr");
        std::cout << a << std::endl;  // do something with a
    }
    xml.set("test.abc.xyz", 12.4543);
    double r = xml.get<double>("test.abc.xyz");  // r is 12.4543
    xml.set("element1.attr", 12, true);  // specify `overwrite = true` to change existing values
    std::cout << xml << std::endl;  // outputs the loaded xml contents to screen
    xml.load("test/testdata/another_file.xml");
    /// [XML usage example]
    std::cout << r << std::endl;  // to not be unused
}

}  // namespace mm
