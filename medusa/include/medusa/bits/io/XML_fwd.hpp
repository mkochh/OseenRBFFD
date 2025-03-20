#ifndef MEDUSA_BITS_IO_XML_FWD_HPP_
#define MEDUSA_BITS_IO_XML_FWD_HPP_

/**
 * @file
 * Declarations for XML I/O utilities.
 *
 * @example test/io/XML_test.cpp
 */

#include <medusa/Config.hpp>
#include <string>
#include <vector>
#include <iosfwd>
#include <rapidxml/rapidxml.hpp>

namespace mm {

/**
 * Class for reading and storing values to XML files. The whole XML file is stored in memory,
 * along with the created document nodes. This class aims to provide simplified support
 * for reading attributes from XML files, usually meant to be used as configurations,
 * similar to `Boost.PropertyTree`.
 * It uses the [RapidXml](http://rapidxml.sourceforge.net) XML library for parsing.
 * For any specifics on parsing XML files, refer to
 * [RapidXml documentation](http://rapidxml.sourceforge.net/manual.html).
 *
 * The intended use is very simple: `.get<int>("element1.element2.attribute")` returns the
 * value of the `attribute`, nested in `element2` inside `element1` inside XML root element.
 * The complementary `.set("element3.attribute", 34.4)` sets the attribute to `34.4`.
 *
 * @warning The attribute path does not contain the root XML element for brevity, as XML standard
 * requires each document to have unique root element and that part of the path would have
 * always been the same.
 *
 * Usage example:
 * @snippet io/XML_test.cpp XML usage example
 * @ingroup io
 */
class XML {
    std::vector<char> file_contents;  ///< Whole XML file contents.
    rapidxml::xml_document<char> doc;  ///< XML document root.

  public:
    /// Creates an XML reader linked to no document.
    XML() = default;
    /**
     * Creates a copy of the XML document.
     * @note This re-parses the whole document. It should be used sparingly, and XML configurations
     * can usually be passed around as `const &` to avoid this.
     */
    XML(const XML&);
    /// Constructs a XML reader by reading the file given by `filename`.
    explicit XML(const std::string& filename);
    /**
     * Loads XML document from a file given by `filename`.
     * Any previous data held in this object instance is lost.
     */
    void load(const std::string& filename);

  private:
    /**
     * Function for opening files, called by XML::XML(const std::string&) and XML::load().
     * The file is opened in this function and closed after this function finishes.
     * @throws Assertion fails if anything is wrong when accessing or reading the file.
     */
    void loadFileHelper(const std::string& file);
    /**
     * Loads the XML document from a stored null-terminated `file_contents`. This modifies the
     * `file_contents` buffer, so the function is not idempotent.
     * @throws rapidxml::parse_error Exception is thrown if the string stored in `file_contents`
     * cannot be parsed to a valid XML document.
     */
    void loadFromStoredContents();
    /// Reads the contents of the attribute specified by `path` as a string.
    char* getString(const std::vector<std::string>& path, const std::string& attribute_name) const;
    /// Writes the contents of string `content` to the attribute specified by `path`.
    void setString(const std::vector<std::string>& path, const std::string& attribute_name,
                   const std::string& content);
    /**
     * Fills the `all_attr` array with all pairs `(path, value)` pairs that are
     * descendants of node `node` at path `path`.
     * @param[in] node Current node, whose descendants will be searched.
     * @param[in] path Path of the given `node`.
     * @param[out] all_attr Container to be filled with `(path, value)` pairs.
     */
    void getAllRecursive(const rapidxml::xml_node<>* node, std::string path,
                         std::vector<std::pair<std::string, std::string>>& all_attr) const;

    /**
     * Splits dot separated path into elements path and attribute name.
     * @return Pair `{{sequence_of_elements}, attribute_name}` representing the path.
     */
    static std::pair<std::vector<std::string>, std::string> splitPath(std::string path);

    /// Join a path into dot separated string.
    static std::string join(const std::vector<std::string>& path);

  public:
    /// Returns `true` if the an attribute specified by `path` exists and `false` otherwise.
    bool exists(const std::string& path) const;

    /**
     * Reads a values from an attribute specified by `path`.
     * @tparam T Type to which the attribute value should be converted. This is done using the
     * stream extraction operator `>>` for type `T`.
     * @param path Dot separated path to the attribute, excluding the root XML tag.
     * @return Value of the attribute converted to type `T`.
     * @throws Assertion fails if the attribute does not exist.
     */
    template <typename T>
    T get(const std::string& path) const;

    /**
     * Saves a value to the attribute pointed to by `path`. If the attribute or path does not exist
     * it is created.
     * @param path Dot separated path to the attribute, excluding root XML tag.
     * It the path does not exist, it is created on the fly.
     * @param value Value to save. It is converted to string before saving using the stream
     * insertion operator `<<` for type `T`.
     * @param overwrite If `true`, possible existing value is overwritten,
     * otherwise an assertion is triggered. Defaults to `false`.
     * @warning If overwrite is false, assertion fails if the attribute exists.
     * This is the default behaviour.
     */
    template <typename T>
    void set(const std::string& path, const T& value, bool overwrite = false);

    /// Returns all pairs `(path, attribute_value)`.
    std::vector<std::pair<std::string, std::string>> getAll() const;

    /// Access to underlying XML root element from RapidXml library.
    rapidxml::xml_document<char>& documentRoot();
    /// Const version of XML::documentRoot.
    const rapidxml::xml_document<char>& documentRoot() const;

    /// Prints the contents of currently loaded XML document.
    friend std::ostream& operator<<(std::ostream& os, const XML& xml);
};


}  // namespace mm

#endif  // MEDUSA_BITS_IO_XML_FWD_HPP_
