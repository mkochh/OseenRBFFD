#include <medusa/bits/utils/memutils.hpp>

#include "gtest/gtest.h"

namespace mm {

//! [deep_copy_unique_ptr def]
class A {
    int a;
  public:
    A(int a) : a(a) {}
    virtual ~A() = default;

    virtual int f() { return a; }
    virtual A* clone() const { return new A(*this); }
};

class B : public A {
    int b;
  public:
    B(int a, int b) : A(a), b(b) {}

    int f() override { return b; }
    B* clone() const override { return new B(*this); }
};
//! [deep_copy_unique_ptr def]

TEST(Utils, MemDeep_copy_unique_ptr) {
    //! [deep_copy_unique_ptr usage]
    deep_copy_unique_ptr<A> pa(new A(3));
    deep_copy_unique_ptr<A> pb(new B(1, 2));

    EXPECT_EQ(3, pa->f());
    EXPECT_EQ(2, pb->f());

    deep_copy_unique_ptr<A> pc = pb;
    EXPECT_EQ(2, pc->f());
    pb.reset(nullptr);
    EXPECT_EQ(2, pc->f());
    pa = deep_copy_unique_ptr<A>(*pc);
    EXPECT_EQ(2, pa->f());
    pa = deep_copy_unique_ptr<A>(static_cast<A>(*pc));
    EXPECT_EQ(1, pa->f());
    pb = pc;
    EXPECT_EQ(2, pb->f());
    pb = B(7, 5);
    EXPECT_EQ(5, pb->f());
    pb = std::move(pc);
    EXPECT_EQ(2, pb->f());
    deep_copy_unique_ptr<A> pd(std::move(pb));
    EXPECT_EQ(2, pd->f());
    //! [deep_copy_unique_ptr usage]
}

TEST(Utils, MemMem2str) {
    EXPECT_EQ("16 B", mem2str(16));
    EXPECT_EQ("1.6 kB", mem2str(1600));
    EXPECT_EQ("1.7 kB", mem2str(1682));
    EXPECT_EQ("1 MB", mem2str(1024 * 1024));
    EXPECT_EQ("5.9 GB", mem2str(static_cast<size_t>(5.5 * 1024 * 1024 * 1024)));
    EXPECT_EQ("5.6 GB", mem2str(static_cast<size_t>(5.6 * 1000 * 1000 * 1000)));
}

}  // namespace mm
