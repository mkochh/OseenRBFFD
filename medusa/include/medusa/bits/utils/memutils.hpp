#ifndef MEDUSA_BITS_UTILS_MEMUTILS_HPP_
#define MEDUSA_BITS_UTILS_MEMUTILS_HPP_

/**
 * @file
 * Declaration of memory related utilities.
 *
 * @example test/utils/memutils_test.cpp
 */

#include <memory>
#include <medusa/Config.hpp>
#include <medusa/bits/utils/assert.hpp>

namespace mm {

/**
 * Unique pointer with polymorphic deep copy semantics. It is copy constructible and
 * copy assignable. Copy construction creates a new pointer, containing a clone
 * of the object, the original pointer pointed to.
 * @tparam T Value type. Must be polymorphically cloneable, i.e.\ have a virtual `clone` method
 * as example below.
 *
 * Usage example:
 * @snippet utils/memutils_test.cpp deep_copy_unique_ptr def
 * @snippet utils/memutils_test.cpp deep_copy_unique_ptr usage
 * @ingroup utils
 */
template <typename T>
class deep_copy_unique_ptr : public std::unique_ptr<T, std::default_delete<T>> {
  public:
    using std::unique_ptr<T, std::default_delete<T>>::unique_ptr;
    using std::unique_ptr<T, std::default_delete<T>>::operator=;
    using std::unique_ptr<T, std::default_delete<T>>::reset;

    /// Construct by polymorphically cloning a given value.
    explicit deep_copy_unique_ptr(const T& v) :
            std::unique_ptr<T, std::default_delete<T>>(v.clone()) {}

    /// Copy by cloning the value of `o`.
    deep_copy_unique_ptr(const deep_copy_unique_ptr& o) :
            std::unique_ptr<T, std::default_delete<T>>() {
        reset(o ? o->clone() : nullptr);
    }

    /// Move construct as `unique_ptr`.
    deep_copy_unique_ptr(deep_copy_unique_ptr&& o) noexcept :
            std::unique_ptr<T, std::default_delete<T>>(std::move(o)) { }

    /// Copy assign by cloning a given value.
    deep_copy_unique_ptr& operator=(const T& v) {
        reset(v.clone());
        return *this;
    }

    /// Copy assign by cloning the value of `o`.
    deep_copy_unique_ptr& operator=(const deep_copy_unique_ptr& o) {
        reset(o ? o->clone() : nullptr);
        return *this;
    }

    /// Move assign as `unique_ptr`
    deep_copy_unique_ptr& operator=(deep_copy_unique_ptr&& o) noexcept {
        std::unique_ptr<T, std::default_delete<T>>::operator=(std::move(o));
        return *this;
    }
};

/**
 * Simple function to help format memory amounts for printing. Takes in number of bytes
 * and returns a human readable representation.
 * @ingroup utils
 */
std::string mem2str(std::size_t bytes);

/**
 * Returns number of bytes the container uses in memory. The container must support `size()`.
 * This does not count the memory that may be allocated by objects stored in the container.
 * Also STL containers like vector may actually have more memory allocated than their size.
 * @ingroup utils
 */
template<typename container_t>
std::size_t mem_used(const container_t& v) {
    return sizeof(v[0]) * v.size();
}

}  // namespace mm

#endif  // MEDUSA_BITS_UTILS_MEMUTILS_HPP_
