#ifndef MEDUSA_BITS_IO_HDF_HPP_
#define MEDUSA_BITS_IO_HDF_HPP_

/**
 * @file
 * Implementation of HDF I/O utilities
 */

#include "HDF_fwd.hpp"
#include <medusa/bits/utils/assert.hpp>
#include <tuple>

namespace mm {

template <class T>
T HDF::readAttribute(const std::string& attr_name) const {
    assert_msg(H5Iis_valid(group), "Group id %d invalid. Did you open a group before reading?",
               group);
    hid_t attr = H5Aopen(group, attr_name.c_str(), H5P_DEFAULT);
    assert_msg(attr >= 0, "Attribute '%s' could not be accessed in group '%s' in file '%s'.",
               attr_name, group_name_, filename_);

    hid_t type = H5Aget_type(attr);
    assert_msg(type >= 0, "Failed getting type of attribute '%s' in group '%s' in file '%s'.",
               attr_name, group_name_, filename_);

    T result;
    herr_t status = H5Aread(attr, type, &result);
    assert_msg(status >= 0, "Failed reading attribute '%s' from group '%s' in file '%s'.",
               attr_name, groupName(), filename_);

    H5Tclose(type);
    H5Aclose(attr);
    return result;
}

template <class T>
void HDF::writeAttribute(const std::string& attr_name, const T& value, const hid_t& type,
                         bool overwrite) const {
    assert_msg(H5Iis_valid(group), "Group id %d invalid. Did you open a group before writing?",
               group);
    if (H5Aexists(group, attr_name.c_str())) {
        if (!overwrite) {
            assert_msg(false, "Attribute '%s' in group '%s' in file '%s' already exists. To "
                              "overwrite its contents use parameter overwrite=true.",
                       attr_name, group_name_, filename_);
            return;
        }
        herr_t status = H5Adelete(group, attr_name.c_str());
        assert_msg(status >= 0, "Failed deleting existing attribute '%s' in group '%s' in "
                                "file '%s' before writing a new one.",
                   attr_name, group_name_, filename_);
    }

    hid_t space = H5Screate(H5S_SCALAR);
    hid_t attr = H5Acreate(group, attr_name.c_str(), type, space, H5P_DEFAULT, H5P_DEFAULT);
    assert_msg(attr >= 0, "Failed creating attribute '%s' in group '%s' in file '%s'.",
               attr_name, group_name_, filename_);
    herr_t status = H5Awrite(attr, type, &value);
    assert_msg(status >= 0, "Failed writing attribute '%s' to group '%s' in file '%s'.",
               attr_name, group_name_, filename_);
    H5Sclose(space);
    H5Aclose(attr);
}

template <typename T>
std::pair<std::vector<hsize_t>, std::vector<T>>
HDF::readLinearArray(const std::string& dataset_name) const {
    assert_msg(H5Iis_valid(group), "Group id %d invalid. Did you open a group before reading?",
               group);
    hid_t dataset = H5Dopen(group, dataset_name.c_str(), H5P_DEFAULT);
    assert_msg(dataset >= 0, "Dataset '%s' could not be accessed in group '%s' in file '%s'.",
               dataset_name, group_name_, filename_);

    hid_t dataspace = H5Dget_space(dataset);
    const int ndims = H5Sget_simple_extent_ndims(dataspace);
    std::vector<hsize_t> dims(ndims);
    H5Sget_simple_extent_dims(dataspace, dims.data(), nullptr);  // read dimension into dims
    hsize_t size = 1;
    for (int d = 0; d < ndims; ++d) size *= dims[d];

    std::vector<T> linear_value(size);
    hid_t type = H5Dget_type(dataset);
    herr_t status = H5Dread(dataset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                            linear_value.data());
    assert_msg(status >= 0, "Failed reading dataset '%s' from group '%s' in file '%s'.",
               dataset_name, group_name_, filename_);
    H5Tclose(type);
    H5Sclose(dataspace);
    H5Dclose(dataset);
    return {dims, linear_value};
}

template <typename T>
std::vector<T> HDF::read1DArray(const std::string& dataset_name) const {
    std::vector<hsize_t> dims;
    std::vector<T> value;
    std::tie(dims, value) = readLinearArray<T>(dataset_name);
    assert_msg(dims.size() == 1, "This function is for one dimensional arrays only, got %d-D "
                                 "array with sizes %s.", dims.size(), dims);
    return value;
}

template <typename T>
std::vector<std::vector<T>> HDF::read2DArray(const std::string& dataset_name) const {
    std::vector<hsize_t> dims;
    std::vector<T> linear_value;
    std::tie(dims, linear_value) = readLinearArray<T>(dataset_name);
    assert_msg(dims.size() == 2, "This function is for two dimensional arrays only, got %d-D "
                                 "array with sizes %s.", dims.size(), dims);
    hsize_t cols = dims[0];
    hsize_t rows = dims[1];

    std::vector<std::vector<T>> value(rows, std::vector<T>(cols));
    for (hsize_t i = 0; i < rows; i++)
        for (hsize_t j = 0; j < cols; j++)
            value[i][j] = linear_value[j * rows + i];

    return value;
}

template <typename T>
std::vector<std::vector<std::vector<T>>> HDF::read3DArray(const std::string& dataset_name) const {
    std::vector<hsize_t> dims;
    std::vector<T> linear_value;
    std::tie(dims, linear_value) = readLinearArray<T>(dataset_name);
    assert_msg(dims.size() == 3, "This function is for three dimensional arrays only, got %d-D "
                                 "array with sizes %s.", dims.size(), dims);
    hsize_t tubes = dims[0];
    hsize_t cols = dims[1];
    hsize_t rows = dims[2];

    std::vector<std::vector<std::vector<T>>> value(rows, std::vector<std::vector<T>>(
            cols, std::vector<T>(tubes)));
    for (hsize_t i = 0; i < rows; i++)
        for (hsize_t j = 0; j < cols; j++)
            for (hsize_t k = 0; k < tubes; k++)
                value[i][j][k] = linear_value[k * cols * rows + j * rows + i];

    return value;
}

template <int dim, typename T>
void HDF::writeLinearArray(const std::string& dataset_name, const T* value,
                           const std::array<hsize_t, dim>& sizes, hid_t type,
                           bool overwrite) const {
    assert_msg(H5Iis_valid(group), "Group id %d invalid. Did you open a group before writing?",
               group);
    hid_t dataset;
    if (H5Lexists(group, dataset_name.c_str(), H5P_DEFAULT)) {
        if (!overwrite) {
            assert_msg(false, "Dataset '%s' in group '%s' in file '%s' already exists. To "
                              "overwrite its contents use parameter overwrite=true.",
                       dataset_name, group_name_, filename_);
            return;
        }
        dataset = H5Dopen(group, dataset_name.c_str(), H5P_DEFAULT);
        assert_msg(dataset >= 0, "Failed opening dataset '%s' in group '%s' in file '%s'.",
                   dataset_name, group_name_, filename_);
        hid_t dataspace = H5Dget_space(dataset);
        const int ndims = H5Sget_simple_extent_ndims(dataspace);
        assert_msg(ndims == dim, "Expected %d-dim array, got %d-dim.", dim, ndims);
        hsize_t cur_dims[dim];  // NOLINT(*)
        H5Sget_simple_extent_dims(dataspace, cur_dims, nullptr);  // read dimension into sizes
        for (int d = 0; d < dim; ++d) {
            assert_msg(cur_dims[d] == sizes[d],
                       "Only data of same old_size can be overwritten, but new dataset has "
                       "old_size '%d' and existing has old_size '%d'.", sizes[d], cur_dims[d]);
        }
        H5Sclose(dataspace);
    } else {
        hid_t dataspace;
        const int rank = dim;
        const hsize_t* max_dims = nullptr;  // same as sizes
        dataspace = H5Screate_simple(rank, &sizes[0], max_dims);
        dataset = H5Dcreate(group, dataset_name.c_str(), type, dataspace,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        assert_msg(dataset >= 0, "Failed creating dataset '%s' in group '%s' in file '%s'.",
                   dataset_name, group_name_, filename_);
        H5Sclose(dataspace);
    }
    int size = 1; for (int d = 0; d < dim; ++d) { size *= sizes[d]; }
    herr_t status = H5Dwrite(dataset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, value);
    assert_msg(status >= 0, "Failed writing dataset '%s' to group '%s' in file '%s'.",
               dataset_name, group_name_, filename_);
    H5Dclose(dataset);
}

template <typename T, typename array_t>
void HDF::write1DArray(const std::string& dataset_name, const array_t& value, hid_t type,
                       bool overwrite) const {
    hsize_t size = value.size();
    std::vector<T> cast_value(size);
    for (hsize_t j = 0; j < size; j++) { cast_value[j] = static_cast<T>(value[j]); }
    writeLinearArray<1>(dataset_name, cast_value.data(), {size}, type, overwrite);
}

template <typename T, class array_t>
void HDF::write2DArray(const std::string& dataset_name, const array_t& value, hid_t type,
                       bool overwrite, std::false_type) const {
    size_t rows = value.size();
    size_t cols = (rows == 0) ? 0 : value[0].size();
    std::vector<T> linear_value(rows * cols);
    for (size_t i = 0; i < rows; i++) {
        assert_msg(static_cast<size_t>(value[i].size()) == cols,
                   "Not all subarrays are of the same size, subarray 0 has size %d, but subarray "
                   "%d has size %d. Maybe use pad to fix raggedness?", cols, i, value[i].size());
        for (size_t j = 0; j < cols; j++) {
            linear_value[j * rows + i] = static_cast<T>(value[i][j]);
        }
    }
    writeLinearArray<2>(dataset_name, linear_value.data(), {cols, rows}, type, overwrite);
}

template <typename T, class array_t>
void HDF::write3DArray(const std::string& dataset_name, const array_t& value, hid_t type,
                       bool overwrite) const {
    size_t rows = value.size();
    size_t cols = (rows == 0) ? 0 : value[0].size();
    size_t tubes = (cols == 0) ? 0 : value[0][0].size();
    std::vector<T> linear_value(rows * cols * tubes);
    for (size_t i = 0; i < rows; i++) {
        assert_msg(static_cast<size_t>(value[i].size()) == cols,
                   "Not all subarrays are of the same size, subarray 0 has size %d, but subarray "
                   "%d has size %d. Maybe use pad to fix raggedness?", cols, i, value[i].size());
        for (size_t j = 0; j < cols; j++) {
            assert_msg(static_cast<size_t>(value[i][j].size()) == tubes,
                       "Not all subarrays are of the same size, subarray (0, 0) has size %d, "
                       "but subarray (%d, %d) has size %d. Maybe use pad to fix raggedness?",
                       tubes, i, j, value[i][j].size());
            for (size_t k = 0; k < tubes; ++k) {
                linear_value[k * cols * rows + j * rows + i] = static_cast<T>(value[i][j][k]);
            }
        }
    }
    writeLinearArray<3>(dataset_name, linear_value.data(), {tubes, cols, rows}, type, overwrite);
}

template <typename array_t>
void HDF::writeIntArray(const std::string& dataset_name, const array_t& value,
                        bool overwrite) const {
    write1DArray<int>(dataset_name, value, H5T_NATIVE_INT, overwrite);
}

template <typename array_t>
void HDF::writeDoubleArray(const std::string& dataset_name, const array_t& value,
                           bool overwrite) const {
    write1DArray<double>(dataset_name, value, H5T_NATIVE_DOUBLE, overwrite);
}

template <typename array_t>
void HDF::writeFloatArray(const std::string& dataset_name, const array_t& value,
                          bool overwrite) const {
    write1DArray<float>(dataset_name, value, H5T_NATIVE_FLOAT, overwrite);
}

template <typename array_t>
void HDF::writeInt2DArray(const std::string& dataset_name, const array_t& value,
                          bool overwrite) const {
    write2DArray<int>(dataset_name, value, H5T_NATIVE_INT, overwrite);
}

template <typename array_t>
void HDF::writeDouble2DArray(const std::string& dataset_name, const array_t& value,
                             bool overwrite) const {
    write2DArray<double>(dataset_name, value, H5T_NATIVE_DOUBLE, overwrite);
}

template <typename array_t>
void HDF::writeFloat2DArray(const std::string& dataset_name, const array_t& value,
                            bool overwrite) const {
    write2DArray<float>(dataset_name, value, H5T_NATIVE_FLOAT, overwrite);
}

template <typename array_t>
void HDF::writeInt3DArray(const std::string& dataset_name, const array_t& value,
                          bool overwrite) const {
    write3DArray<int>(dataset_name, value, H5T_NATIVE_INT, overwrite);
}

template <typename array_t>
void HDF::writeDouble3DArray(const std::string& dataset_name, const array_t& value,
                             bool overwrite) const {
    write3DArray<double>(dataset_name, value, H5T_NATIVE_DOUBLE, overwrite);
}

template <typename array_t>
void HDF::writeFloat3DArray(const std::string& dataset_name, const array_t& value,
                            bool overwrite) const {
    write3DArray<float>(dataset_name, value, H5T_NATIVE_FLOAT, overwrite);
}

template <typename conf_t>
void HDF::writeXML(const std::string& name, const conf_t& conf, bool overwrite) {
    std::string group = group_name_;
    assert_msg(!group.empty(), "Open a group before writing.");
    if (group.back() == '/') openGroup(group+name);
    else openGroup(group+'/'+name);
    std::vector<std::pair<std::string, std::string>> data = conf.getAll();
    for (const auto& kv : data) {
        try {
            std::string::size_type num_read;
            double x = std::stod(kv.second, &num_read);
            if (num_read < kv.second.size()) {  // some characters were not read, it's a string
                throw std::invalid_argument(kv.second);
            }
            writeDoubleAttribute(kv.first, x, overwrite);
        } catch (const std::invalid_argument&) {
            writeStringAttribute(kv.first, kv.second, overwrite);
        }
    }
    openGroup(group);
}

template <typename SparseMatrixType>
void HDF::writeSparseMatrix(const std::string& name, SparseMatrixType& matrix, bool one_based,
                            bool overwrite) {
    std::vector<std::array<double, 3>> triplets(matrix.nonZeros());
    int c = 0;
    for (int k = 0; k < matrix.outerSize(); ++k) {
        for (typename SparseMatrixType::InnerIterator it(matrix, k); it; ++it) {
            triplets[c][0] = one_based+it.row();
            triplets[c][1] = one_based+it.col();
            triplets[c][2] = it.value();
            ++c;
        }
    }
    writeDouble2DArray(name, triplets, overwrite);
}

template <typename domain_t>
void HDF::writeDomain(const std::string& name, const domain_t& domain, bool overwrite) {
    std::string group = group_name_;
    assert_msg(!group.empty(), "Open a group before writing.");
    if (group.back() == '/') openGroup(group+name);
    else openGroup(group+'/'+name);
    writeDouble2DArray("pos", domain.positions(), overwrite);
    writeIntAttribute("N", domain.size(), overwrite);
    writeIntArray("types", domain.types(), overwrite);
    writeIntArray("bmap", domain.bmap(), overwrite);
    writeDouble2DArray("normals", domain.normals(), overwrite);
    openGroup(group);
}

template <typename timer_t>
void HDF::writeTimer(const std::string& name, const timer_t& timer, bool overwrite) {
    std::string group = group_name_;
    assert_msg(!group.empty(), "Open a group before writing.");
    if (group.back() == '/') openGroup(group+name);
    else openGroup(group+'/'+name);
    std::vector<std::string> labels = timer.labels();
    int size = labels.size();
    if (size == 0) return;
    for (int i = 1; i < size; ++i) {
        writeDoubleAttribute(labels[i-1]+"-"+labels[i], timer.duration(labels[i-1], labels[i]),
                             overwrite);
    }
    writeDoubleAttribute("total", timer.duration(labels.front(), labels.back()), overwrite);
    openGroup(group);
}

}  // namespace mm

#endif  // MEDUSA_BITS_IO_HDF_HPP_
