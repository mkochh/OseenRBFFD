#ifndef MEDUSA_BITS_IO_HDF_FWD_HPP_
#define MEDUSA_BITS_IO_HDF_FWD_HPP_

/**
 * @file
 * Declarations of CSV I/O utilities
 *
 * @example test/io/HDF_test.cpp
 */

#include <medusa/Config.hpp>
#include <medusa/bits/utils/print.hpp>
#include <hdf5.h>
#include <string>
#include <vector>
#include <fstream>
#include <ostream>
#include <type_traits>

/// @cond
namespace Eigen {
template <typename Derived>
class MatrixBase;

template<typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
class Matrix;
}

template <template <typename...> class base, typename derived>
struct is_base_of_template_impl {
    template<typename... Ts>
    static constexpr std::true_type  test(const base<Ts...> *);
    static constexpr std::false_type test(...);
    using type = decltype(test(std::declval<derived*>()));
};

template <template <typename...> class base, typename derived>
using is_base_of_template = typename is_base_of_template_impl<base, derived>::type;
/// @endcond

namespace mm {

/**
 * Simplified [HDF5](https://portal.hdfgroup.org/display/support) I/O utilities.
 * This class is a wrapper around the
 * [HDF5 C library](https://portal.hdfgroup.org/pages/viewpage.action?pageId=50073943),
 * for simplified reading and writing.
 *
 * This class represents a HDF reader and writer: it knows which file (open or closed)
 * and which group (open or closed) are the ones that are to be read from/written to.
 * Any calls to `read*` or `write*` methods write to that group in that file.
 * The file and group can be opened, closed or changed during the lifetime of the object,
 * with user's responsibility being that they are open when writing to/reading from it.
 *
 * In addition atomic() method can be used to perform `open/read/close` or `open/write/close`
 * cycle in one line.
 *
 * Raw C identifiers are exposed with getFileID() and getGroupID() methods and can be used to
 * perform any action not directly supported by this class.
 *
 * All `write` methods create the datasets or attributes required.
 * They also have an `overwrite` parameter which is set to `false` by default.
 * This means that any attempt to write to an existing object will fail an assertion to prevent
 * accidental deletions of existing data.
 * This can be disabled by setting `overwrite = true`, allowing one to overwrite existing data.
 *
 * All open objects are closed upon destruction.
 *
 * @note To enable support for Eigen types, HDF_Eigen.hpp must be included.
 *
 * Usage example:
 * @snippet io/HDF_test.cpp HDF usage example
 * @ingroup io
 */
class HDF {
    std::string filename_;  ///< Current file name.
    std::string group_name_;  ///< Current group name.

    hid_t file;  ///< Currently open file identifier.
    hid_t group;  ///< Currently open group identifier.

  public:
    /// Possible file opening modes.
    enum Mode : unsigned {
        APPEND = 256,  ///< Appends to the existing contents, if any exist.
        DESTROY = 128,  ///< Removes old contents, if any exist.
        READONLY = 64   ///< Read only open, the file must exist.
    };
  private:
    /// Convert HDF::Mode enum to string.
    static std::string str(Mode mode);

  public:
    /// Construct an empty HDF reader.
    HDF() : filename_(), group_name_(), file(-1), group(-1) {}
    /**
     * Construct a HDF reader and open `filename` with given `mode`.
     * It also opens group `/` by default, so it is immediately ready for reading / writing.
     * @throws Assertion fails if mode was `READONLY` and the file can not be read.
     */
    explicit HDF(std::string filename, Mode mode = APPEND) :
            filename_(std::move(filename)), group_name_("/"), file(-1), group(-1) {
        file = openFileHelper(filename_, mode);
        group = H5Gopen(file, "/", H5P_DEFAULT);
    }
    /// Closes all opened objects, by using close(). @sa close
    ~HDF() { close(); }

  private:
    /// Opens given file with given mode. @sa open, openFile, openGroup
    static hid_t openFileHelper(const std::string& filename, Mode mode);

  public:
    /**
     * Opens given file and group.
     * @param filename Name of the file to open.
     * @param group_name Name of the group to open.
     * @param mode Specifies the mode to open the file in, see HDF::Mode.
     */
    void open(const std::string& filename, const std::string& group_name = "/",
              Mode mode = APPEND) {
        openFile(filename, mode);
        openGroup(group_name);
    }

    /**
     * Opens a given file with a given mode. Any previously opened file is closed.
     * No group is opened.
     */
    void openFile(const std::string& filename, Mode mode = APPEND) {
        closeFile();
        filename_ = filename;
        file = openFileHelper(filename, mode);
    }

    /// Reopens the current file.
    void reopenFile(Mode mode = APPEND) { openFile(filename_, mode); }

    /**
     * Open given group, closing any previously opened group().
     * @param group_name Name of the group to open. The name must be `/` separated and can
     * be either absolute or relative. Absolute names such as `/a/b/c` open groups
     * staring from `/` group and relative names such as `a/b/c` open groups staring
     * from current group. The groups are created automatically if they do not exist.
     */
    void openGroup(std::string group_name);
    /// Reopens current group.
    void reopenGroup() {
        if (!H5Iis_valid(group)) { openGroup(group_name_); }
    }
    /// Reopens current file and group.
    void reopen(Mode mode = APPEND) { reopenFile(mode); reopenGroup(); }
    /// Returns `true` is currently specified file is open and `false` otherwise.
    bool isFileOpen() const { return H5Iis_valid(file) != 0; }
    /// Returns `true` is currently specified group is open and `false` otherwise.
    bool isGroupOpen() const { return H5Iis_valid(group) != 0; }
    /// Flush opened file, if it is open.
    void flush() const;
    /// Closes current group, if open.
    void closeGroup() const;
    /// Closes current file and all objects associated with it (e.g.\ the group), if open.
    void closeFile() const;
    /// Closes all open objects. Alias of closeFile().
    void close() const { closeFile(); }
    /// Get current filename.
    const std::string& filename() const { return filename_; }
    /// Alias of filename().
    const std::string& fileName() const { return filename(); }
    /// Sets new filename without opening the new file. Any previously opened file is closed.
    void setFilename(const std::string& filename) { closeFile(); filename_ = filename; }
    /// Alias of setFilename().
    void setFileName(const std::string& filename) { setFilename(filename); }
    /// Get current group name.
    const std::string& groupName() const { return group_name_; }
    /// Set new group name without opening the new group. Any previously opened group is closed.
    void setGroupName(const std::string& group_name);

    /**
     * Read attribute given by `attr_name`.
     * @tparam T The type to which to read the attribute.
     * @param attr_name Name of the attribute
     * @return Attribute value, converted to type `T`.
     * @throws Assertion fails if the attribute does not exists or cannot be read.
     * @warning If the attribute type is not compatible with type `T`, the function
     * will still attempt to read it the behavior is undefined. Most often the return value
     * is garbage.
     * @note This function is specialized for `std::string`.
     */
    template <class T>
    T readAttribute(const std::string& attr_name) const;

    /**
     * Writes given attribute with given value and type to file.
     * @param attr_name Name of the attribute.
     * @param value Value to be written
     * @param type HDF5 type of the attribute to be written.
     * @param overwrite If overwriting existing attributes is allowed. If an attribute with the
     * same name is attempted to be written, the existing attribute is deleted and a new one is
     * created instead.
     */
    template <class T>
    void writeAttribute(const std::string& attr_name, const T& value, const hid_t& type,
                        bool overwrite = false) const;

    /**
     * Read given dataset into linear `dim`-D array.
     * @param dataset_name Name of the dataset.
     * @return A vector containing all dataset elements in col-major fashion.
     * @throws Assertion fails if the dataset does not exists or cannot be read.
     * @warning If the element type is not compatible with type `T`, the function
     * will still attempt to read it the behavior is undefined. Most often the return value
     * is garbage.
     */
    template <typename T>
    std::pair<std::vector<hsize_t>, std::vector<T>>
    readLinearArray(const std::string& dataset_name) const;

    /**
     * Read 1D dataset given by `dataset_name`.
     * @param dataset_name Name of the dataset.
     * @return A vector containing all dataset elements.
     * @throws Assertion fails if the dataset does not exists or cannot be read.
     * @warning If the element type is not compatible with type `T`, the function
     * will still attempt to read it the behavior is undefined. Most often the return value
     * is garbage.
     */
    template<typename T>
    std::vector<T> read1DArray(const std::string& dataset_name) const;

    /**
     * Read given dataset into 2D array of vectors.
     * @param dataset_name Name of the dataset.
     * @return  A vector of vectors containing all dataset elements.
     * @throws Assertion fails if the dataset does not exists or cannot be read.
     * @warning If the element type is not compatible with type `T`, the function
     * will still attempt to read it the behavior is undefined. Most often the return value
     * is garbage.
     */
    template<typename T>
    std::vector<std::vector<T>> read2DArray(const std::string& dataset_name) const;

    /**
     * Read given dataset into 3D array of vectors.
     * @param dataset_name Name of the dataset.
     * @return  A vector of vectors containing all dataset elements.
     * @throws Assertion fails if the dataset does not exists or cannot be read.
     * @warning If the element type is not compatible with type `T`, the function
     * will still attempt to read it the behavior is undefined. Most often the return value
     * is garbage.
     */
    template<typename T>
    std::vector<std::vector<std::vector<T>>> read3DArray(const std::string& dataset_name) const;

    /**
     * Writes a `dim`-dimensional dataset stored in 1D array in col-major fashion.
     * @tparam dim Dimension (HDF rank) of the dataset, 2 for matrices, 3 for tensors.
     * @param dataset_name Name of the dataset.
     * @param value Pointer to values to be written, stored in a linear col-major fashion.
     * The number of values should be equal to the product of `sizes`.
     * @param sizes Sizes of the dataset along each dimension.
     * @param type HDF5 type of the dataset.
     * @param overwrite See write2DArray()
     */
    template<int dim, typename T>
    void writeLinearArray(const std::string& dataset_name, const T* value,
                          const std::array<hsize_t, dim>& sizes, hid_t type, bool overwrite) const;

    /**
     * Writes given 1D array with given value and type to file.
     * @tparam T Type to which the elements of value are converted to.
     * @tparam array_t Type of the given array. It must support `operator[]` and `.size()`.
     * @param dataset_name Name of the attribute.
     * @param value Value to be written
     * @param type HDF5 type of the elements to be written.
     * @param overwrite If overwriting existing datasets is allowed. If an attribute with the
     * same name is attempted to be written, it must have the same dimensions and type
     * as current value. Otherwise the overwrite will fail.
     */
    template<typename T, typename array_t>
    void write1DArray(const std::string& dataset_name, const array_t& value,
                      hid_t type, bool overwrite) const;

    /**
     * Writes given 2D array of given value and type to file.
     * @tparam T Type of `array_t` elements
     * @tparam array_t Type of the given 2D array. It must support `.size()`, `[i].size()` and
     * `[i][j]` operations.
     * @param dataset_name Name of the dataset.
     * @param value Value to be written.
     * @param type HDF5 type of the elements to be written.
     * @param overwrite If overwriting existing datasets is allowed. If an attribute with the
     * same name is attempted to be written, it must have the same dimensions and type
     * as current value. Otherwise the overwrite will fail.
     */
    template<typename T, class array_t>
    void write2DArray(const std::string& dataset_name, const array_t& value, hid_t type,
                      bool overwrite) const {
        write2DArray<T>(dataset_name, value, type, overwrite,
                        is_base_of_template<Eigen::MatrixBase, array_t>());
    }

  private:
    /// Overload for static_assert for Eigen due to different `size` method.
    template<typename T, class array_t>
    void write2DArray(const std::string&, const array_t&, hid_t, bool, std::false_type) const;

    /// Overload for static_assert for Eigen due to different `size` method.
    template<typename T, class array_t>
    void write2DArray(const std::string&, const array_t&, hid_t, bool, std::true_type) const {
        static_assert(!std::is_same<T, T>::value,
                "This function produces wrong results for Eigen matrices. Use writeEigen instead.");
    }

  public:
    /**
     * Writes given 3D array of given value and type to file.
     * @tparam T Type of `array_t` elements
     * @tparam array_t Type of the given 2D array. It must support `.size()`, `[i].size()`,
     * `[i][j]`, `[i][j].size()` and `[i][j][k]` operations.
     * @param dataset_name Name of the dataset.
     * @param value Value to be written.
     * @param type HDF5 type of the elements to be written.
     * @param overwrite If overwriting existing datasets is allowed. If an attribute with the
     * same name is attempted to be written, it must have the same dimensions and type
     * as current value. Otherwise the overwrite will fail.
     */
    template<typename T, class array_t>
    void write3DArray(const std::string& dataset_name, const array_t& value, hid_t type,
                      bool overwrite) const;

    /// Reads `int` attribute. @sa readAttribute
    int readIntAttribute(const std::string& attr_name) const;
    /// Reads `bool` attribute. @sa readAttribute
    bool readBoolAttribute(const std::string& attr_name) const;
    /// Reads `double` attribute. @sa readAttribute
    double readDoubleAttribute(const std::string& attr_name) const;
    /// Reads `float` attribute. @sa readAttribute
    float readFloatAttribute(const std::string& attr_name) const;
    /**
     * Reads non-null-terminated `string` attribute. Any null characters are read as part of the
     * string.
     * @sa readAttribute
     */
    std::string readStringAttribute(const std::string& attr_name) const;
    /// Reads null-terminated `string` attribute. @sa readAttribute
    std::string readNullTerminatedStringAttribute(const std::string& attr_name) const;

    /// Write `int` attribute. @sa writeAttribute
    void writeIntAttribute(const std::string& attr_name, int value, bool overwrite = false) const;
    /// Write `bool` attribute. @sa writeAttribute
    void writeBoolAttribute(const std::string& attr_name, bool value, bool overwrite = false) const;
    /// Write `double` attribute. @sa writeAttribute
    void writeDoubleAttribute(const std::string& attr_name, double value,
                              bool overwrite = false) const;
    /// Write `float` attribute. @sa writeAttribute
    void writeFloatAttribute(const std::string& attr_name, float value,
                             bool overwrite = false) const;
    /**
     * Write `string` attribute. The string is written as an non-null-terminater array of bytes.
     * It may contain null characters in the middle.
     * @sa writeAttribute
     */
    void writeStringAttribute(const std::string& attr_name, const std::string& value,
                              bool overwrite = false) const;

    /// Read given dataset into vector of `int`s. @sa readArray
    std::vector<int> readIntArray(const std::string& dataset_name) const;
    /// Read given dataset into vector of `double`s. @sa readArray
    std::vector<double> readDoubleArray(const std::string& dataset_name) const;
    /// Read given dataset into vector of `float`s. @sa readArray
    std::vector<float> readFloatArray(const std::string& dataset_name) const;

    /// Write given value as a 1D array of `int`s. @sa write1DArray
    template <typename array_t>
    void writeIntArray(const std::string& dataset_name, const array_t& value,
                       bool overwrite = false) const;
    /// Write given value as a 1D array of `double`s. @sa write1DArray
    template <typename array_t>
    void writeDoubleArray(const std::string& dataset_name, const array_t& value,
                          bool overwrite = false) const;
    /// Write given value as a 1D array of `floats`s. @sa write1DArray
    template <typename array_t>
    void writeFloatArray(const std::string& dataset_name, const array_t& value,
                         bool overwrite = false) const;

    /// Read given dataset into 2D vector of `int`s. @sa read2DArray
    std::vector<std::vector<int>> readInt2DArray(const std::string& dataset_name) const;
    /// Read given dataset into 2D vector of `double`s. @sa read2DArray
    std::vector<std::vector<double>> readDouble2DArray(const std::string& dataset_name) const;
    /// Read given dataset into 2D vector of `float`s. @sa read2DArray
    std::vector<std::vector<float>> readFloat2DArray(const std::string& dataset_name) const;

    /// Write given value as a 2D array of `int`s. @sa write2DArray
    template <typename array_t>
    void writeInt2DArray(const std::string& dataset_name, const array_t& value,
                         bool overwrite = false) const;
    /// Write given value as a 2D array of `double`s. @sa write2DArray
    template <typename array_t>
    void writeDouble2DArray(const std::string& dataset_name, const array_t& value,
                            bool overwrite = false) const;
    /// Write given value as a 2D array of `float`s. @sa write2DArray
    template <typename array_t>
    void writeFloat2DArray(const std::string& dataset_name, const array_t& value,
                           bool overwrite = false) const;

    /// Read given dataset into 3D vector of `int`s. @sa read2DArray
    std::vector<std::vector<std::vector<int>>> readInt3DArray(
            const std::string& dataset_name) const;
    /// Read given dataset into 3D vector of `double`s. @sa read2DArray
    std::vector<std::vector<std::vector<double>>> readDouble3DArray(
            const std::string& dataset_name) const;
    /// Read given dataset into 3D vector of `float`s. @sa read2DArray
    std::vector<std::vector<std::vector<float>>> readFloat3DArray(
            const std::string& dataset_name) const;

    /// Write given value as a 3D array of `int`s. @sa write3DArray
    template <typename array_t>
    void writeInt3DArray(const std::string& dataset_name, const array_t& value,
                         bool overwrite = false) const;
    /// Write given value as a 3D array of `double`s. @sa write3DArray
    template <typename array_t>
    void writeDouble3DArray(const std::string& dataset_name, const array_t& value,
                            bool overwrite = false) const;
    /// Write given value as a 3D array of `float`s. @sa write3DArray
    template <typename array_t>
    void writeFloat3DArray(const std::string& dataset_name, const array_t& value,
                           bool overwrite = false) const;

    /**
     * Reads data as an Eigen matrix of type `scalar_t`.
     * @tparam scalar_t The of the data to read. Set to `double` by default.
     * @param dataset_name Name of the dataset to read from.
     * @warning If the supplied data type does not match the actual type of the data,
     * the behavior is undefined. Most likely, the returned matrix will have garbage contents.
     * @note To use this function, the file HDF_Eigen.hpp must be included.
     */
    template <typename scalar_t = double>
    Eigen::Matrix<scalar_t, -1, -1, 0, -1, -1> readEigen(const std::string& dataset_name) const;

    /**
     * Write given Eigen expression as a 2D array of its type (either `int`, `double` or `float`).
     * @note To use this function, the file HDF_Eigen.hpp must be included.
     */
    template <typename Derived>
    void writeEigen(const std::string& dataset_name, const Eigen::MatrixBase<Derived>& value,
                    bool overwrite = false) const;

    /**
     * Writes a sparse matrix.
     * Matrix `matrix` is written as a `N x 3` matrix, where the columns
     * represent vectors `I`, `J`, `V`, such that `M(I(i), J(i)) = V(i)`.
     * Matrix values are stored as doubles of native size.
     * @warning By default, the indexes of the matrix elements are saved in 1-based format.
     * @param name Attribute name.
     * @param matrix A sparse matrix.
     * @param one_based Whether to store the indices as one or zero based. One based indexes are
     * ready to be read by Matlab's `spconvert`.
     * @param overwrite See write2DArray().
     */
    template <typename SparseMatrixType>
    void writeSparseMatrix(const std::string& name, SparseMatrixType& matrix, bool one_based = true,
                           bool overwrite = false);

    /**
     * Writes given domain discretization to file. A group with given name is created
     * and 4 datasets and 1 attribute are stored inside:
     * - `N`: int attribute storing number of nodes in the domain,
     * - `pos`: `N x dim` double array of positions,
     * - `types`: int array of length `N` representing domain types,
     * - `bmap`: int array of length `N` mapping nodes to boundary nodes
     *    (see DomainDiscretization::bmap),
     * - `normals`: double array of domain normals for boundary nodes.
     * After the function exists, the same group is opened as before the call.
     *
     * @param name Name of the group in which the domain is stored.
     * @param domain Domain to write to file.
     * @param overwrite If `true`, allows overwriting of existing entries.
     * @sa DomainDiscretization
     */
    template <typename domain_t>
    void writeDomain(const std::string& name, const domain_t& domain, bool overwrite = false);

    /**
     * Writes given timer to file. A group with given name is created
     * and durations (in seconds) between consecutive labels are writen as
     * double attributes named `label1-label2`. Additionally, total time
     * between first and last label (in seconds) is stored in the double
     * attribute named `total`. After the function exists, the same group is opened as before the
     * call. If timer has 0 labels, only the group is created. If timer has one label,
     * attribute `total` with value `0.0` is written.
     * @param name Name of the group in which the timer is stored.
     * @param timer Timer to write to file.
     * @param overwrite If `true`, allows overwriting of existing entries.
     * @sa Timer
     */
    template <typename timer_t>
    void writeTimer(const std::string& name, const timer_t& timer, bool overwrite = false);

    /**
     * Writes given XML object to file. All attributes in the file, given by XML::getAll,
     * are written in their dot-separated full-path format. They are written either as
     * doubles (if whole value can be parsed as a `double`) or strings otherwise.
     * After the function exists, the same group is opened as before the call.
     * @param name Name of the group in which the XML attributes are stored.
     * @param conf XML object.
     * @param overwrite If `true`, allows overwriting of existing entries.
     * @sa XML
     */
    template <typename conf_t>
    void writeXML(const std::string& name, const conf_t& conf, bool overwrite = false);

    /**
     * Allows for "atomic" read and write operations to HDF5 files.
     *
     * This means that read and write functions open the file and group,
     * before reading/writing and close them afterwards. The filename and group name must be set
     * before use.
     * @sa setFilename, setGroupName
     */
    HDF atomic() const;

    /// Return raw HDF C file identifier.
    hid_t getFileID() const { return file; }
    /// Return raw HDF C group identifier.
    hid_t getGroupID() const { return group; }

  private:
    /// Callback required by HDF `Literate` function when getting a list of members.
    static herr_t memberIterateCallback(hid_t loc_id, const char* name,
                                        const H5L_info_t*, void* data);

  public:
    /// Holds categorized names for all members of a group.
    struct Members {
        std::vector<std::string> groups,  ///< Names of subgroups.
                datasets,  ///< Names of datasets.
                datatypes,  ///< Names of datatypes.
                unknowns;  ///< Names of unknown objects.

        /// Output this type in a user friendly way for quick inspection.
        friend std::ostream& operator<<(std::ostream& os, const Members& members) {
            os << "groups: " << members.groups << "\ndatasets: " << members.datasets
               << "\ndatatypes: " << members.datatypes << "\nunknowns: " << members.unknowns;
            return os;
        }
    };

    /**
     * Returns an object with categorized names of all members of current group. The names within
     * each category are guaranteed to be in alphabetical order.
     */
    Members members() const;

    /// Returns names of all subgroups in current group in alphabetical order. @sa members
    std::vector<std::string> groups() const { return members().groups; }

    /// Returns names of all datasets in current group in alphabetical order. @sa members
    std::vector<std::string> datasets() const { return members().datasets; }

    /// Output basic info about HDF reader, such as current group and file.
    friend std::ostream& operator<<(std::ostream& os, const HDF& hdf) {
        os << "HDF reader for file '" << hdf.filename()
           << "' with current group '" << hdf.groupName() << "'. ";
        os << "The file is " << (hdf.isFileOpen() ? "open" : "closed")
           << " and the group is " << (hdf.isGroupOpen() ? "open" : "closed") << '.';
        return os;
    }
};

}  // namespace mm

#endif  // MEDUSA_BITS_IO_HDF_FWD_HPP_
