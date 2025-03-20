#ifndef MEDUSA_BITS_TYPES_RANGE_HPP_
#define MEDUSA_BITS_TYPES_RANGE_HPP_

/**
 * @file
 * Implementation of Range.
 */

#include "Range_fwd.hpp"
#include <algorithm>

namespace mm {

template<class T, class Allocator>
Range<T, Allocator>& Range<T, Allocator>::operator=(std::initializer_list<value_type> lst) {
    this->resize(lst.size());
    int i = 0;
    for (const_reference x : lst) {
        operator[](i++) = x;
    }
    return *this;
}

template<class T, class Allocator>
Range<T, Allocator>& Range<T, Allocator>::operator=(const value_type& x) {
    for (size_type i = 0; i < size(); ++i) {
        operator[](i) = x;
    }
    return *this;
}

template<class T, class Allocator>
typename Range<T, Allocator>::reference Range<T, Allocator>::operator[](size_type i) {
    assert_msg(0 <= i && i < size(), "Index %d out of range [%d, %d) when accessing Range "
            "for write.", i, 0, size());
    return std::vector<T>::operator[](i);
}

template<class T, class Allocator>
typename Range<T, Allocator>::const_reference Range<T, Allocator>::operator[](size_type i) const {
    assert_msg(0 <= i && i < size(), "Index %d out of range [%d, %d) when accessing Range "
            "for read/write.", i, 0, size());
    return std::vector<T>::operator[](i);
}

template<class T, class Allocator>
typename Range<T, Allocator>::RangeView Range<T, Allocator>::operator[](const indexes_t& indexes) {
    for (size_type idx : indexes) {
        assert_msg(0 <= idx && idx < size(), "Index %d out of range [%d, %d) when using "
                "multiindexed read-write access.", idx, 0, size());
    }
    return {*this, indexes};
}

template<class T, class Allocator>
typename Range<T, Allocator>::ConstRangeView
Range<T, Allocator>::operator[](const indexes_t& indexes) const {
    for (size_type idx : indexes) {
        assert_msg(0 <= idx && idx < size(), "Index %d out of range [%d, %d) when using "
                "multiindexed read access.", idx, 0, size());
    }
    return {*this, indexes};
}

template<class T, class Allocator>
Range<T, Allocator>& Range<T, Allocator>::append(const Range<T, Allocator>& rng) {
    this->insert(this->end(), rng.begin(), rng.end());
    return *this;
}

template<class T, class Allocator>
Range<T, Allocator> Range<T, Allocator>::join(const Range<T, Allocator>& rng) const {
    Range ret = *this;
    ret.insert(ret.end(), rng.begin(), rng.end());
    return ret;
}

template<class T, class Allocator>
void Range<T, Allocator>::remove(indexes_t indexes) {
    std::sort(indexes.begin(), indexes.end());
    indexes.resize(std::unique(indexes.begin(), indexes.end()) - indexes.begin());
    int sz = indexes.size();
    auto it = std::vector<T, Allocator>::begin();
    auto cur = std::vector<T, Allocator>::begin();
    int to_remove = 0;
    int c = 0;
    while (cur != std::vector<T, Allocator>::end()) {
        if (to_remove < sz && c == indexes[to_remove]) {
            ++to_remove;
        } else {
            *it++ = std::move(*cur);
        }
        ++c;
        ++cur;
    }
    std::vector<T, Allocator>::resize(it - std::vector<T, Allocator>::begin());
}

template<class T, class Allocator>
template<class Predicate>
indexes_t Range<T, Allocator>::filter(const Predicate& predicate) const {
    indexes_t ret;
    for (size_type i = 0; i < size(); ++i)
        if (predicate(operator[](i))) ret.push_back(i);
    return ret;
}

template<class T, class Allocator>
template<typename UnaryOp>
auto Range<T, Allocator>::map(UnaryOp op) -> Range<decltype(op(this->operator[](0)))> const {
    Range<decltype(op(this->operator[](0)))> r;
    r.reserve(size());
    std::transform(Range<T>::begin(), Range<T>::end(), std::back_inserter(r), op);
    return r;
}

template<class T, class Allocator>
indexes_t Range<T, Allocator>::operator<(const value_type& v) const {
    return filter([&](const value_type& t) { return t < v; });
}

template<class T, class Allocator>
indexes_t Range<T, Allocator>::operator>(const value_type& v) const {
    return filter([&](const value_type& t) { return t > v; });
}

template<class T, class Allocator>
indexes_t Range<T, Allocator>::operator<=(const value_type& v) const {
    return filter([&](const value_type& t) { return t <= v; });
}

template<class T, class Allocator>
indexes_t Range<T, Allocator>::operator>=(const value_type& v) const {
    return filter([&](const value_type& t) { return t >= v; });
}

template<class T, class Allocator>
indexes_t Range<T, Allocator>::operator==(const value_type& v) const {
    return filter([&](const value_type& t) { return t == v; });
}

template<class T, class Allocator>
indexes_t Range<T, Allocator>::operator!=(const value_type& v) const {
    return filter([&](const value_type& t) { return t != v; });
}

/// @cond
template <class T, class Allocator>
template <typename V>
Range<T, Allocator> Range<T, Allocator>::seq(V n) {
    return Range<T, Allocator>::seq(static_cast<V>(0), n);
}

template <class T, class Allocator>
template <typename V>
Range<T, Allocator> Range<T, Allocator>::seq(V start, V stop) {
    return Range<T, Allocator>::seq(start, stop, static_cast<V>(1));
}

template <class T, class Allocator>
template <typename V, typename D>
Range<T, Allocator> Range<T, Allocator>::seq(V start, V stop, D step) {
    assert_msg(step != 0, "Step must be nonzero.");
    assert_msg((stop - start) / step >= 0, "Step %d must be in the direction from start %d to "
                                           "stop %d.", step, start, stop);
    Range<T> ret; ret.reserve((stop - start) / step);
    for (V i = start; (step > 0) ? i < stop : i > stop; i += step) {
        ret.push_back(i);
    }
    return ret;
}
/// @endcond

}  // namespace mm

#endif   // MEDUSA_BITS_TYPES_RANGE_HPP_
