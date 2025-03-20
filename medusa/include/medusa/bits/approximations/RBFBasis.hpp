#ifndef MEDUSA_BITS_APPROXIMATIONS_RBFBASIS_HPP_
#define MEDUSA_BITS_APPROXIMATIONS_RBFBASIS_HPP_

/**
 * @file
 * Implementation of RBF basis.
 */

#include <medusa/bits/utils/assert.hpp>
#include "RBFBasis_fwd.hpp"

namespace mm {

/// @cond
template <class vec_t, class RBFType>
typename  RBFBasis<vec_t, RBFType>::scalar_t RBFBasis<vec_t, RBFType>::eval(
        int index, const vector_t& point, const std::vector<vector_t>& support) const {
    assert_msg(0 <= index && index < size_, "Basis index %d out of range [0, %d)", index, size_);
    assert_msg(static_cast<int>(support.size()) >= size_,
               "Not enough support points given for this basis size, got %d support points for "
               "basis size %d.", support.size(), size_);
    return rbf_((point-support[index]).squaredNorm());
}

template <class RBFType, class vec_t>
template <typename operator_t>
typename vec_t::scalar_t RBFBasis<RBFType, vec_t>::evalOp(int index, const vector_t& point,
        operator_t op, const std::vector<vector_t>& support, scalar_t scale) const {
    return op.apply(*this, index, point, support, scale);
}
/// @endcond

template <class RBFType, class vec_t>
typename vec_t::scalar_t RBFBasis<RBFType, vec_t>::evalOp(
        int index, const vector_t& point, Lap<dim> op, const std::vector<vector_t>& support,
        scalar_t scale) const {
    assert_msg(0 <= index && index < size_, "Basis index %d out of range [0, %d)", index, size_);
    assert_msg(static_cast<int>(support.size()) >= size_,
               "Not enough support points given for this basis size, got %d support points for "
               "basis size %d.", support.size(), size_);
    return rbf_((point - support[index]).squaredNorm(), op)/scale/scale;
}

template <class RBFType, class vec_t>
typename vec_t::scalar_t RBFBasis<RBFType, vec_t>::evalOp(
        int index, const vector_t& point, Der1<dim> op, const std::vector<vector_t>& support,
        scalar_t scale) const {
    assert_msg(0 <= index && index < size_, "Basis index %d out of range [0, %d)", index, size_);
    assert_msg(static_cast<int>(support.size()) >= size_,
               "Not enough support points given for this basis size, got %d support points for "
               "basis size %d.", support.size(), size_);
    vec_t p = point - support[index];
    return 2*p[op.var]*rbf_(p.squaredNorm(), 1)/scale;
}

template <class RBFType, class vec_t>
typename vec_t::scalar_t RBFBasis<RBFType, vec_t>::evalOp(
        int index, const vector_t& point, Der2<dim> op, const std::vector<vector_t>& support,
        scalar_t scale) const {
    assert_msg(0 <= index && index < size_, "Basis index %d out of range [0, %d)", index, size_);
    assert_msg(static_cast<int>(support.size()) >= size_,
               "Not enough support points given for this basis size, got %d support points for "
               "basis size %d.", support.size(), size_);
    vec_t p = point - support[index];
    scalar_t r2 = p.squaredNorm();
    double res = 4*p[op.var1]*p[op.var2]*rbf_(r2, 2);
    if (op.var1 == op.var2) {
        res += 2*rbf_(r2, 1);
    }
    return res/scale/scale;
}

template <class vec_t, class RBFType>
typename RBFBasis<vec_t, RBFType>::scalar_t RBFBasis<vec_t, RBFType>::evalAt0(
        int index, const std::vector<vector_t>& support) const {
    assert_msg(0 <= index && index < size_, "Basis index %d out of range [0, %d)", index, size_);
    assert_msg(static_cast<int>(support.size()) >= size_,
               "Not enough support points given for this basis size, got %d support points for "
               "basis size %d.", support.size(), size_);
    return rbf_(support[index].squaredNorm());
}

template <class RBFType, class vec_t>
typename vec_t::scalar_t
RBFBasis<RBFType, vec_t>::evalOpAt0(int index, const Lap<dim>& lap,
        const std::vector<vector_t>& support, scalar_t scale) const {
    assert_msg(0 <= index && index < size_, "Basis index %d out of range [0, %d)", index, size_);
    assert_msg(static_cast<int>(support.size()) >= size_,
               "Not enough support points given for this basis size, got %d support points for "
               "basis size %d.", support.size(), size_);
    return rbf_(support[index].squaredNorm(), lap)/scale/scale;
}

template <class RBFType, class vec_t>
typename vec_t::scalar_t RBFBasis<RBFType, vec_t>::evalOpAt0(int index, const Der1<dim>& der1,
        const std::vector<vector_t>& support, scalar_t scale) const {
    assert_msg(0 <= index && index < size_, "Basis index %d out of range [0, %d)", index, size_);
    assert_msg(static_cast<int>(support.size()) >= size_,
               "Not enough support points given for this basis size, got %d support points for "
               "basis size %d.", support.size(), size_);
    assert_msg(der1.var >= 0, "Index of derived variable %d should be higher or equal to 0.",
               der1.var);
    assert_msg(der1.var < dim, "Index of derived variable %d should be lower"
                               " than the number of dimensions %d", der1.var, dim);

    const auto& s = support[index];
    scalar_t r2 = s.squaredNorm();
    scalar_t v = rbf_(r2, 1);
    v *= -2*s[der1.var];
    return v/scale;
}

template <class RBFType, class vec_t>
typename vec_t::scalar_t RBFBasis<RBFType, vec_t>::evalOpAt0(int index, const Der2<dim>& der2,
        const std::vector<vector_t>& support, scalar_t scale) const {
    assert_msg(0 <= index && index < size_, "Basis index %d out of range [0, %d)", index, size_);
    assert_msg(static_cast<int>(support.size()) >= size_,
               "Not enough support points given for this basis size, got %d support points for "
               "basis size %d.", support.size(), size_);
    assert_msg(der2.var1 >= 0 && der2.var2 >=0, "Indexes of derived variables %d and %d should"
                                                " be both higher or equal to 0.",
               der2.var1, der2.var2);
    assert_msg(der2.var1 < dim && der2.var2 < dim , "Index of derived variable %d  and %d should"
              " both be lower than the number of dimensions %d", der2.var1, der2.var2, dim);

    auto& s = support[index];
    scalar_t r2 = s.squaredNorm();
    int d1 = der2.var1;
    int d2 = der2.var2;
    scalar_t v = 4*s[d1]*s[d2]*rbf_(r2, 2);
    if (d1 == d2) {
        v += 2*rbf_(r2, 1);
    }
    return v/scale/scale;
}

/// Output basic info about given basis.
template <typename V, typename R>
std::ostream& operator<<(std::ostream& os, const RBFBasis<V, R>& m) {
    return os << "RBF basis " << m.dim << "D, " << m.size() << "x " << m.rbf();
}

}  // namespace mm

#endif  // MEDUSA_BITS_APPROXIMATIONS_RBFBASIS_HPP_
