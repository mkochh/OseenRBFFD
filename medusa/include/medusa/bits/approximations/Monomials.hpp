#ifndef MEDUSA_BITS_APPROXIMATIONS_MONOMIALS_HPP_
#define MEDUSA_BITS_APPROXIMATIONS_MONOMIALS_HPP_

/**
 * @file
 * Implementation of Monomial basis.
 */

#include "Monomials_fwd.hpp"
#include <medusa/bits/utils/assert.hpp>
#include <medusa/bits/utils/numutils.hpp>
#include <iostream>
#include "Operators_fwd.hpp"

namespace mm {

template <class vec_t>
Monomials<vec_t>::Monomials(int order)  {
    assert_msg(-1 <= order, "Requested monomials of negative order %d.", order);
    setFromPowers(generatePowers(order, dim));
}

template <class vec_t>
std::vector<std::vector<int>> Monomials<vec_t>::generatePowers(int max_order, int dim) {
    if (max_order < 0) return {};
    if (dim == 0) return {{}};
    std::vector<std::vector<int>> powers;
    for (int i = 0; i <= max_order; ++i) {
        auto other = generatePowers(max_order - i, dim-1);
        for (const auto& p : other) {
            powers.push_back({i});
            for (int o : p) {
                powers.back().push_back(o);
            }
        }
    }
    return powers;
}

template <class vec_t>
Monomials<vec_t> Monomials<vec_t>::tensorBasis(int order) {
    assert_msg(-1 <= order, "Requested monomials of negative order %d.", order);
    if (order == -1) return Monomials<vec_t>();
    std::vector<std::vector<int>> powers;
    Vec<int, dim> counter = 0;
    Vec<int, dim> counts = order+1;
    do {
        powers.push_back(std::vector<int>(counter.begin(), counter.end()));
    } while (incrementCounter(counter, counts));
    return Monomials<vec_t>(powers);
}

template <class vec_t>
void Monomials<vec_t>::setFromPowers(const std::vector<std::vector<int>>& powers) {
    powers_.resize(dim, powers.size());
    int size = powers.size();
    for (int i = 0; i < size; ++i) {
        const std::vector<int>& p = powers[i];
        assert_msg(p.size() == dim, "Monomial size %s does not match dimension %d", p.size(), dim);
        for (int j = 0; j < dim; ++j) {
            powers_(j, i) = p[j];
        }
    }
}

template <class vec_t>
typename vec_t::scalar_t Monomials<vec_t>::eval(
        int index, const vec_t& point, const std::vector<vector_t>& /* support */) const {
    assert_msg(0 <= index && index < size(), "Monomial at index %d does not exist. "
            "Index must be in range [0, %d).", index, size());
    scalar_t result(1.0);
    for (int i = 0; i < dim; i++) {
        result *= ipow(point[i], powers_(i, index));
    }
    return result;
}

template <class vec_t>
typename vec_t::scalar_t Monomials<vec_t>::evalOp(
        int index, const vector_t& point, Lap<dim>, const std::vector<vector_t>&,
        scalar_t scale) const {
    scalar_t result(0.0);
    for (int d = 0; d < dim; ++d) {
        if (powers_(d, index) <= 1) continue;
        scalar_t der2(1.0);
        for (int i = 0; i < d; i++) {
            der2 *= ipow(point[i], powers_(i, index));
        }
        der2 *= powers_(d, index) * (powers_(d, index) - 1) * ipow(point[d], powers_(d, index)-2);
        for (int i = d+1; i < dim; i++) {
            der2 *= ipow(point[i], powers_(i, index));
        }
        result += der2;
    }
    return result/scale/scale;
}

template <class vec_t>
typename vec_t::scalar_t Monomials<vec_t>::evalOp(
        int index, const vector_t& point, Der1<dim> op, const std::vector<vector_t>&,
        scalar_t scale) const {
    if (powers_(op.var, index) == 0) return 0;
    scalar_t result(1.0);
    for (int i = 0; i < op.var; i++) {
        result *= ipow(point[i], powers_(i, index));
    }
    result *= powers_(op.var, index)*ipow(point[op.var], powers_(op.var, index)-1);
    for (int i = op.var+1; i < dim; i++) {
        result *= ipow(point[i], powers_(i, index));
    }
    return result/scale;
}

template <class vec_t>
typename vec_t::scalar_t Monomials<vec_t>::evalOp(
        int index, const vector_t& point, Der2<dim> op, const std::vector<vector_t>&,
        scalar_t scale) const {
    if (op.var1 == op.var2) {
        if (powers_(op.var1, index) <= 1) return 0;
        scalar_t result(1.0);
        for (int i = 0; i < op.var1; i++) {
            result *= ipow(point[i], powers_(i, index));
        }
        result *= powers_(op.var1, index) * (powers_(op.var1, index)-1) *
                ipow(point[op.var1], powers_(op.var1, index)-2);
        for (int i = op.var1+1; i < dim; i++) {
            result *= ipow(point[i], powers_(i, index));
        }
        return result/scale/scale;
    } else {
        if (powers_(op.var1, index) == 0 || powers_(op.var2, index) == 0 ) return 0;
        scalar_t result(1.0);
        for (int i = 0; i < op.var1; i++) {
            result *= ipow(point[i], powers_(i, index));
        }
        result *= powers_(op.var1, index)*ipow(point[op.var1], powers_(op.var1, index)-1);
        for (int i = op.var1+1; i < op.var2; i++) {
            result *= ipow(point[i], powers_(i, index));
        }
        result *= powers_(op.var2, index)*ipow(point[op.var2], powers_(op.var2, index)-1);
        for (int i = op.var2+1; i < dim; i++) {
            result *= ipow(point[i], powers_(i, index));
        }
        return result/scale/scale;
    }
}

template <class vec_t>
typename vec_t::scalar_t Monomials<vec_t>::evalOp(int index, const vector_t& point,
        Derivative<dim> d, const std::vector<vector_t>&, scalar_t scale) const {
    assert_msg(0 <= index && index < size(), "Monomial at index %d does not exist. "
                                             "Index must be in range [0, %d).", index, size());
    int totaldeg = 0;
    for (int x : d.orders) {
        assert_msg(x >= 0, "Derivative of negative order %d requested.", x);
        totaldeg += x;
    }
    scalar_t result(1.0);
    for (int i = 0; i < dim; i++) {
        if (d.orders[i] > powers_(i, index)) {
            result = 0;
            break;
        }
        result *= ipow(point[i], powers_(i, index) - d.orders[i]);
        for (int j = powers_(i, index); j > powers_(i, index) - d.orders[i]; j--) {
            result *= j;
        }
    }
    return result / ipow(scale, totaldeg);
}

template <class vec_t>
typename vec_t::scalar_t Monomials<vec_t>::evalAt0(int index, const std::vector<vector_t>&) const {
    assert_msg(0 <= index && index < size(), "Monomial at index %d does not exist. "
                                             "Index must be in range [0, %d).", index, size());
    return powers_.col(index).sum() == 0;
}

/// @cond
template <class vec_t>
template <typename operator_t>
typename vec_t::scalar_t Monomials<vec_t>::evalOp(int index, const vector_t& point, operator_t op,
        const std::vector<vector_t>& support, scalar_t scale) const {
    return op.apply(*this, index, point, support, scale);
}
/// @endcond

template <class vec_t>
typename vec_t::scalar_t Monomials<vec_t>::evalOpAt0(int index, const Lap<dim>&,
        const std::vector<vector_t>&, scalar_t scale) const {
    assert_msg(0 <= index && index < size(), "Monomial at index %d does not exist. "
                                             "Index must be in range [0, %d).", index, size());
    scalar_t result = 0;
    for (int d = 0; d < vec_t::dim; ++d) {
        if (powers()(d, index) != 2) continue;
        scalar_t r = 2.0;
        for (int i = 0; i < vec_t::dim; ++i) {
            int p = powers()(i, index);
            if (i != d && p != 0) { r = 0; break; }
        }
        result += r;
    }
    return result / (scale*scale);
}

template <class vec_t>
typename vec_t::scalar_t Monomials<vec_t>::evalOpAt0(int index, const Der1<dim>& der1,
        const std::vector<vector_t>& , scalar_t scale) const {
    assert_msg(der1.var >= 0, "Index of derived variable %d should be higher or equal to 0.",
            der1.var);
    assert_msg(der1.var < dim, "Index of derived variable %d should be lower"
                               " than the number of dimensions %d", der1.var, dim);
    assert_msg(0 <= index && index < size(), "Monomial at index %d does not exist. "
                                             "Index must be in range [0, %d).", index, size());
    int d = der1.var;
    scalar_t r = 0;
    if (powers()(d, index) == 1) {
        r = 1;
        for (int i = 0; i < vec_t::dim; ++i) {
            int p = powers()(i, index);
            if (i != d && p != 0) r=0;
        }
    }
    return r / scale;
}

template <class vec_t>
typename vec_t::scalar_t Monomials<vec_t>::evalOpAt0(int index, const Der2<dim>& der2,
        const std::vector<vector_t>& , scalar_t scale) const {
    assert_msg(der2.var1 >= 0 && der2.var2 >=0, "Index of derived variables %d and %d should"
                                                " be both higher or equal to 0.",
               der2.var1, der2.var2);
    assert_msg(der2.var1 < dim && der2.var2 < dim , "Indexes of derived variables %d  and %d should"
              " both be lower than the number of dimensions %d", der2.var1, der2.var2, dim);
    assert_msg(0 <= index && index < size(), "Monomial at index %d does not exist. "
                                             "Index must be in range [0, %d).", index, size());
    int d1 = der2.var1;
    int d2 = der2.var2;
    scalar_t r = 0;

    if (d1 == d2) {
        if (powers()(d1, index) == 2) {
            r = 2;
            for (int i = 0; i < vec_t::dim; ++i) {
                int p = powers()(i, index);
                if (i != d1 && p != 0) r = 0;
            }
        }
    } else {
        if (powers()(d1, index) == 1 && powers()(d2, index) == 1) {
            r = 1;
            for (int i = 0; i < vec_t::dim; ++i) {
                int p = powers()(i, index);
                if (i != d1 && i != d2 && p != 0) r = 0;
            }
        }
    }
    return r / (scale * scale);
}

template <class vec_t>
typename vec_t::scalar_t Monomials<vec_t>::evalOpAt0(int index, const Derivative<dim>& der,
        const std::vector<vector_t>& support, scalar_t scale) const {
    return evalOp(index, 0.0, der, support, scale);
}

/// Output basic info about given Monomial basis.
template <class V>
std::ostream& operator<<(std::ostream& os, const Monomials<V>& m) {
    return os << "Monomials " << m.dim << "D: "
              << m.powers_.transpose() << ", number of functions = " << m.size();
}

}  // namespace mm

#endif  // MEDUSA_BITS_APPROXIMATIONS_MONOMIALS_HPP_
