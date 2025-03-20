#ifndef MEDUSA_BITS_IO_HDF_EIGEN_HPP_
#define MEDUSA_BITS_IO_HDF_EIGEN_HPP_

/**
 * @file
 * Implemetation of HDF IO with Eigen types.
 */

#include <Eigen/Core>
#include "HDF_fwd.hpp"
#include <medusa/bits/utils/assert.hpp>

namespace mm {

/// @cond
template <typename scalar_t>
Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic> HDF::readEigen(
        const std::string& dataset_name) const {
    assert_msg(H5Iis_valid(group), "Group id %d invalid. Did you open a group before reading?",
               group);
    hid_t dataset = H5Dopen(group, dataset_name.c_str(), H5P_DEFAULT);
    assert_msg(dataset >= 0, "Dataset '%s' could not be accessed in group '%s' in file '%s'.",
               dataset_name, group_name_, filename_);

    hid_t dataspace = H5Dget_space(dataset);
    const int ndims = H5Sget_simple_extent_ndims(dataspace);
    assert_msg(ndims <= 2, "This function is for 1 and 2 dimensional arrays only.");
    hsize_t dims[2] = {1, 1};
    H5Sget_simple_extent_dims(dataspace, dims, nullptr);  // read dimension into dims
    hsize_t cols = dims[0], rows = dims[1];
    Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic> M(rows, cols);

    hid_t type = H5Dget_type(dataset);
    herr_t status = H5Dread(dataset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, M.data());
    assert_msg(status >= 0, "Failed reading dataset '%s' from group '%s' in file '%s'.",
               dataset_name, group_name_, filename_);

    H5Tclose(type);
    H5Sclose(dataspace);
    H5Dclose(dataset);
    return M;
}
/// @endcond

template <typename Derived>
void HDF::writeEigen(const std::string& dataset_name, const Eigen::MatrixBase<Derived>& value,
                     bool overwrite) const {
    typedef typename Eigen::MatrixBase<Derived>::Scalar scalar_t;
    if (value.IsRowMajor) {
        writeEigen(dataset_name,
                   Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>(value),
                   overwrite);
        return;
    }
    hsize_t cols = value.cols();
    hsize_t rows = value.rows();
    if (std::is_same<double, scalar_t>::value) {
        writeLinearArray<2>(dataset_name, value.eval().data(), {cols, rows},
                            H5T_NATIVE_DOUBLE, overwrite);
    } else if (std::is_same<float, scalar_t>::value) {
        writeLinearArray<2>(dataset_name, value.eval().data(), {cols, rows},
                            H5T_NATIVE_FLOAT, overwrite);
    } else if (std::is_same<int, scalar_t>::value) {
        writeLinearArray<2>(dataset_name, value.eval().data(), {cols, rows},
                            H5T_NATIVE_INT, overwrite);
    } else {
        assert_msg(false, "Only float, int and double types are supported.");
    }
}

}  // namespace mm

#endif  // MEDUSA_BITS_IO_HDF_EIGEN_HPP_
