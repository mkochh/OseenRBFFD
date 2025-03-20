#ifndef MEDUSA_BITS_TYPES_RANGE_FWD_HPP_
#define MEDUSA_BITS_TYPES_RANGE_FWD_HPP_

/**
 * @file
 * Declaration of Range.
 *
 * @example test/types/Range_test.cpp
 */

#include <medusa/Config.hpp>
#include <medusa/bits/utils/assert.hpp>
#include "traits.hpp"
#include <memory>
#include <vector>
#include <iosfwd>

namespace mm {

/**
 * An extension of `std::vector<T>`  to support additional useful operations.
 * It only adds new methods, no new members, so destructing via base pointer is still safe
 * (even though this usage is not intended).
 *
 * Usage example:
 * @snippet types/Range_test.cpp Range usage example
 * @ingroup types
 */
template<class T, class Allocator = std::allocator<T>>
class Range : public std::vector<T, Allocator> {
  public:
    /// This container's value type
    typedef typename std::vector<T, Allocator>::value_type value_type;
    /// This container's reference to value type
    typedef typename std::vector<T, Allocator>::reference reference;
    /// This container's const reference to value type
    typedef typename std::vector<T, Allocator>::const_reference const_reference;
    /// This container's size type
    typedef int size_type;

    /**
     * This class represents a non contiguous view to a Range, allowing for
     * read and write operations.
     */
    class RangeView {
        Range<T>& receiver;  ///< Reference to object we are viewing / modifying.
        const indexes_t& modifier;  ///< List of indexes of elements to modify.

        /// Constructor.
        RangeView(Range<T>& receiver_, const indexes_t& modifier_)
                : receiver(receiver_), modifier(modifier_) {}

      public:
        RangeView(RangeView&) = delete;  ///< Disallow copying.
        RangeView(RangeView&&) = delete;  ///< Disallow moving.
        RangeView& operator=(RangeView&) = delete;  ///< Disallow copying.
        RangeView& operator=(RangeView&&) = delete;  ///< Disallow moving.

        /// Cast to underlying container type.
        Range<T> asRange() const {
            Range<T> ret(size());
            int sz = modifier.size();
            for (int i = 0; i < sz; ++i) ret[i] = operator[](i);
            return ret;
        }

        /// Multiindex assignment: `a[{1, 2, 3}] = Range<int>({1, 2, 3});`.
        void operator=(const Range<T>& rhs) {
            assert_msg(rhs.size() == size(),
                       "Container sizes must match in multiindexed assignment, but my size is %d "
                       "and assigned size is %d.", size(), rhs.size());
            for (size_type i = 0; i < size(); ++i) operator[](i) = rhs[i];
        }

        /// Multiindex value assignment: `a = 4;`.
        void operator=(const value_type& x) {
            for (size_type i = 0; i < size(); ++i)
                operator[](i) = x;
        }

        /// Write access to sub-container elements.
        reference operator[](size_type i) { return receiver[modifier[i]]; }

        /// Read access to sub-container elements.
        const_reference operator[](size_type i) const { return receiver[modifier[i]]; }

        /// Size of the sub-container.
        size_type size() const { return static_cast<int>(modifier.size()); }

        friend class Range;

        /// Output a RangeView.
        friend std::ostream& operator<<(std::ostream& os, const RangeView& c) {
            return os << c.asRange();
        }
    };

    /**
    * This class represents a non contiguous view to a Range, allowing for
    * read-only operations.
    */
    class ConstRangeView {
        const Range<T>& receiver;  ///< Reference to object we are viewing.
        const indexes_t& modifier;  ///< List of indexes of elements to modify.

        /// Constructor.
        ConstRangeView(const Range<T>& receiver_, const indexes_t& modifier_)
                : receiver(receiver_), modifier(modifier_) {}

      public:
        ConstRangeView(ConstRangeView&) = delete;  ///< Disallow copying.
        ConstRangeView(ConstRangeView&&) = delete;  ///< Disallow moving.
        ConstRangeView& operator=(ConstRangeView&) = delete;  ///< Disallow copying.
        ConstRangeView& operator=(ConstRangeView&&) = delete;  ///< Disallow moving.

        /// Cast to underlying container type.
        Range<T> asRange() const {
            Range<T> ret(size());
            int sz = modifier.size();
            for (int i = 0; i < sz; ++i) ret[i] = operator[](i);
            return ret;
        }

        /// Read access to sub-container elements.
        const_reference operator[](size_type i) const { return receiver[modifier[i]]; }

        /// Size of the sub-container.
        size_type size() const { return static_cast<int>(modifier.size()); }

        friend class Range;

        /// Output a ConstRangeView.
        friend std::ostream& operator<<(std::ostream& os, const ConstRangeView& c) {
            return os << c.asRange();
        }
    };

  public:
    using std::vector<T, Allocator>::vector;
    using std::vector<T, Allocator>::operator=;
    Range(const Range& o) = default;  ///< Default copy constructor.
    Range() : std::vector<T, Allocator>() {}  ///< Default constructor.
    Range(Range&& o) noexcept = default;  ///< Default move constructor.
    Range& operator=(const Range& o) = default;  ///< Default copy assignment.
    Range& operator=(Range&& o) noexcept = default;  ///< Default move assignment.
    /// Copy construct from `std::vector`.
    Range(const std::vector<T, Allocator>& o) : std::vector<T, Allocator>(o) {}
    /// Move Construct from `std::vector`.
    Range(std::vector<T, Allocator>&& o) noexcept : std::vector<T, Allocator>(std::move(o)) {}
    /// Construct from `RangeView`.
    Range(const Range::RangeView& o) : Range(o.asRange()) {}
    /// Construct from `ConstRangeView`.
    Range(const Range::ConstRangeView& o) : Range(o.asRange()) {}

    /// Returns range `{0, ..., n-1}`.
    template <typename V>
    static Range seq(V n);
    /// Returns range with integers in interval `[start, stop)`.
    template <typename V>
    static Range seq(V start, V stop);
    /// Returns range with integers in interval `[start, stop)` with step `step`.
    template <typename V, typename D>
    static Range seq(V start, V stop, D step);

    /// Assign from initializer list: `a = {1, 2, 3};`.
    Range& operator=(std::initializer_list<value_type> lst);

    /// Assign a single value to all elements: `a = 4;`.
    Range& operator=(const value_type& x);

    /// Overload vector's `[]` to assert parameter.
    reference operator[](size_type i);

    /// Overload vector's `[]` to assert parameter.
    const_reference operator[](size_type i) const;

    /// Multi-indexed access for writing.
    RangeView operator[](const indexes_t& indexes);

    /// Multi-indexed access for reading.
    ConstRangeView operator[](const indexes_t& indexes) const;

    /// Returns number of elements.
    /// @note This function returns signed type, contrary to `std::vector::size`.
    int size() const { return static_cast<int>(std::vector<T, Allocator>::size()); }

    /// Append all elements of `rng` to self.
    Range& append(const Range& rng);
    /// Operator version of Range::append.
    Range& operator+=(const Range& rng) { return append(rng); }

    /// Return new copy containing this range's elements followed by all elements of `rng`.
    Range join(const Range& rng) const;

    /// Remove elements wih given indexes.
    void remove(indexes_t indexes);

    /**
     * Returns list of indexes for which predicate returns `true`.
     * @tparam Predicate Any callable object, e.g.\ lambda func, functional,
     * class with `operator()` defined.
     * Example:
     * @code
     * auto idxs = a.filter([](double v){ return 2.3 < v && v < 6.4; })];
     * @endcode
     * returns all elements of a that are between 2.3 and 6.4.
     */
    template<class Predicate>
    indexes_t filter(const Predicate& predicate) const;

    /**
     * Returns a new range, obtained by applying `op` to all elements of this Range.
     * @tparam UnaryOp Any callable object.
     * @param op Operation mapping `T -> Out` type.
     * @return `[a, b, c].map(f) -> [f(a), f(b), f(c)]`
     */
    template<typename UnaryOp>
    auto map(UnaryOp op) -> Range<decltype(op(this->operator[](0)))> const;

    /// Returns list of indexes of elements that are less than `v`.
    indexes_t operator<(const value_type& v) const;

    /// Returns list of indexes of elements that are greater than `v`.
    indexes_t operator>(const value_type& v) const;

    /// Returns list of indexes of elements that are less or equal to `v`.
    indexes_t operator<=(const value_type& v) const;

    /// Returns list of indexes of elements that are greater or equal to `v`.
    indexes_t operator>=(const value_type& v) const;

    /// Returns list of indexes of elements that are equal to `v`.
    indexes_t operator==(const value_type& v) const;

    /// Returns list of indexes of elements that are not equal to `v`.
    indexes_t operator!=(const value_type& v) const;
};

/// The scalar_type trait definition for Range.
template <typename T> struct scalar_type<Range<T>> {
    typedef typename Range<T>::value_type type;  ///< Underlying scalar type.
};
/// The vector_type trait definition for Range.
template <typename T> struct vector_type<Range<T>> {
    typedef typename Range<T>::value_type type;  ///< Underlying scalar type.
};
/// The scalar_type trait definition for std::vector.
template <typename T> struct scalar_type<std::vector<T>> {
    typedef typename std::vector<T>::value_type type;  ///< Underlying vector type.
};
/// The vector_type trait definition for std::vector.
template <typename T> struct vector_type<std::vector<T>> {
    typedef typename std::vector<T>::value_type type;  ///< Underlying vector type.
};

/**
 * Concatenate two vectors.
 * @code
 * std::vector<int> a = {1, 2, 3};
 * std::vector<int> b = {4, 5, 6};
 * std::vector<int> c = a + b;  // c = {1, 2, 3, 4, 5, 6};
 * @endcode
 */
template<class T, class Allocator = std::allocator<T>>
std::vector<T, Allocator> operator+(const std::vector<T, Allocator>& v1,
                                    const std::vector<T, Allocator>& v2) {
    std::vector<T, Allocator> ret = v1;
    ret.insert(ret.end(), v2.begin(), v2.end());
    return ret;
}

}  // namespace mm

#endif  // MEDUSA_BITS_TYPES_RANGE_FWD_HPP_
