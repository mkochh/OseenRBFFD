#ifndef MEDUSA_BITS_IO_XML_HPP_
#define MEDUSA_BITS_IO_XML_HPP_

#include "XML_fwd.hpp"
#include <medusa/bits/utils/assert.hpp>
#include <sstream>
#include <iomanip>

/**
 * @file
 * Implementation of XML I/O utilities.
 */

namespace mm {

template <typename T>
T XML::get(const std::string& path) const {
    std::vector<std::string> path_elements;
    std::string attribute_name;
    std::tie(path_elements, attribute_name) = splitPath(path);
    char* value = getString(path_elements, attribute_name);
    std::stringstream ss(value);
    T result;
    ss >> result;
    return result;
}

/// String specialization.
template <>
std::string XML::get<std::string>(const std::string& path) const;

template <typename T>
void XML::set(const std::string& path, const T& value, bool overwrite) {
    std::vector<std::string> path_elements;
    std::string attribute_name;
    std::tie(path_elements, attribute_name) = splitPath(path);
    if (!overwrite) {
        assert_msg(!exists(path), "Attribute on path '%s' already exists with value '%s'. "
                                  "Set overwrite=true to overwrite its value.",
                   path, getString(path_elements, attribute_name));
    }
    std::stringstream ss;
    ss << std::setprecision(16) << value;
    setString(path_elements, attribute_name, ss.str());
}

}  // namespace mm

#endif  // MEDUSA_BITS_IO_XML_HPP_
