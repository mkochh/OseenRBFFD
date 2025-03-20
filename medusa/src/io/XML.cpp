#include <medusa/bits/io/XML.hpp>
#include <rapidxml/rapidxml.hpp>
#include <rapidxml/rapidxml_print.hpp>
#include <medusa/bits/utils/assert.hpp>
#include <fstream>
#include <iostream>
#include <cstring>

/**
 * @file
 * Implementation of XML I/O utilities.
 */

namespace mm {

void XML::loadFileHelper(const std::string& file) {
    std::ifstream ifs(file);
    assert_msg(ifs.good(), "Failed opening file '%s' with error: %s.", file, strerror(errno));
    ifs.seekg(0, std::ios::end);
    std::size_t length = ifs.tellg();
    ifs.seekg(0, std::ios::beg);
    file_contents.resize(length+1, '\0');
    ifs.read(file_contents.data(), length);
    ifs.close();
    loadFromStoredContents();
}

void XML::loadFromStoredContents() {
    try {
        doc.parse<0>(file_contents.data());
    } catch (const rapidxml::parse_error& e) {
        std::cerr << "Error parsing file contents '"
                  << static_cast<const char*>(file_contents.data()) << "'." << std::endl;
        throw e;
    }
}

XML::XML(const std::string& filename) { loadFileHelper(filename); }
XML::XML(const XML& xml) {
    rapidxml::print(std::back_inserter(file_contents), xml.doc);
    file_contents.push_back('\0');
    loadFromStoredContents();
}
void XML::load(const std::string& filename) { loadFileHelper(filename); }
rapidxml::xml_document<char>& XML::documentRoot() { return doc; }
const rapidxml::xml_document<char>& XML::documentRoot() const { return doc; }

/// Prints the contents of currently loaded XML document.
std::ostream& operator<<(std::ostream& os, const XML& xml) {
    os << "XML document: '\n";
    rapidxml::print(os, xml.doc, 0);
    return os << '\'';
}

std::pair<std::vector<std::string>, std::string> XML::splitPath(std::string path) {
    assert_msg(!path.empty(), "Path name should not be empty.");
    std::size_t last_dot_idx = path.find_last_of('.');
    if (last_dot_idx == std::string::npos) {
        return {{}, path};
    }
    std::string attribute_name = path.substr(last_dot_idx+1);
    assert_msg(!attribute_name.empty(), "Path '%s' contains an empty attribute name.", path);

    path = path.substr(0, last_dot_idx);
    std::vector<std::string> split_path;
    size_t idx = 0;
    do {
        std::size_t new_idx = path.find('.', idx);
        if (new_idx == std::string::npos) { new_idx = path.size(); }
        assert_msg(new_idx - idx > 0, "Path '%s' contains an empty element name.", path);
        split_path.push_back(path.substr(idx, new_idx-idx));
        idx = new_idx + 1;
    } while (idx <= path.size());
    return {split_path, attribute_name};
}

std::string XML::join(const std::vector<std::string>& path) {
    if (path.empty()) { return ""; }
    std::string r = path[0];
    for (std::string::size_type i = 1; i < path.size(); ++i) {
        r += '.';
        r += path[i];
    }
    return r;
}

char* XML::getString(const std::vector<std::string>& path,
                     const std::string& attribute_name) const {
    rapidxml::xml_node<char>* root = doc.first_node();
    for (const auto& name : path) {
        root = root->first_node(name.c_str(), name.size());
        assert_msg(root != nullptr, "Element '%s' from path '%s' not found.",
                   name, join(path));
    }
    auto* attr = root->first_attribute(attribute_name.c_str(), attribute_name.size());
    assert_msg(attr != nullptr, "Attribute '%s', specified by path '%s' not found.",
               attribute_name, join(path));
    return attr->value();
}

void XML::setString(const std::vector<std::string>& path, const std::string& attribute_name,
                    const std::string& content) {
    rapidxml::xml_node<char>* root = doc.first_node();
    for (const auto& name : path) {
        auto* child = root->first_node(name.c_str(), name.size());
        if (child == nullptr) {
            char* node_name = doc.allocate_string(name.c_str());
            child = doc.allocate_node(rapidxml::node_element, node_name);
            root->append_node(child);
        }
        root = child;
    }
    auto* attr = root->first_attribute(attribute_name.c_str(), attribute_name.size());
    if (attr == nullptr) {
        char* attr_name = doc.allocate_string(attribute_name.c_str());
        attr = doc.allocate_attribute(attr_name);
        root->append_attribute(attr);
    }
    assert_msg(attr != nullptr, "Attribute '%s' not found on path '%s'.",
               attribute_name, join(path));
    char* attr_content = doc.allocate_string(content.c_str());
    attr->value(attr_content);
}

/// @cond
template <>
std::string XML::get<std::string>(const std::string& path) const {
    std::vector<std::string> path_elements;
    std::string attribute_name;
    std::tie(path_elements, attribute_name) = splitPath(path);
    return getString(path_elements, attribute_name);
}

template <>
void XML::set<std::string>(const std::string& path, const std::string& value, bool overwrite) {
    std::vector<std::string> path_elements;
    std::string attribute_name;
    std::tie(path_elements, attribute_name) = splitPath(path);
    if (!overwrite) {
        assert_msg(!exists(path), "Attribute on path '%s' already exists with value '%s'. "
                                  "Set overwrite=true to overwrite its value.",
                   path, getString(path_elements, attribute_name));
    }
    setString(path_elements, attribute_name, value);
}
/// @endcond

std::vector<std::pair<std::string, std::string>> XML::getAll() const {
    std::vector<std::pair<std::string, std::string>> all_attributes;
    getAllRecursive(doc.first_node(), "", all_attributes);
    return all_attributes;
}

void XML::getAllRecursive(const rapidxml::xml_node<>* node, std::string path,
                          std::vector<std::pair<std::string, std::string>>& all_attr) const {
    for (const auto* a = node->first_attribute(); a; a = a->next_attribute()) {
        all_attr.emplace_back(path+a->name(), a->value());
    }
    for (const rapidxml::xml_node<>* n = node->first_node(); n; n = n->next_sibling()) {
        getAllRecursive(n, path + n->name() + ".", all_attr);
    }
}

bool XML::exists(const std::string& path) const {
    rapidxml::xml_node<char>* root = doc.first_node();
    std::vector<std::string> path_elements;
    std::string attribute_name;
    std::tie(path_elements, attribute_name) = splitPath(path);
    for (const auto& element_name : path_elements) {
        root = root->first_node(element_name.c_str(), element_name.size());
        if (!root) return false;
    }
    auto* attr = root->first_attribute(attribute_name.c_str(), attribute_name.size());
    return attr != nullptr;
}

}  // namespace mm
