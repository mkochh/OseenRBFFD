#include <medusa/bits/io/HDF.hpp>
#include <medusa/bits/io/HDF_fwd.hpp>


/**
 * @file
 * Implementation of HDF I/O utilities
 */

namespace mm {

std::string HDF::str(HDF::Mode mode)  {
    switch (mode) {
        case APPEND: return "APPEND";
        case DESTROY: return "DESTROY";
        case READONLY: return "READONLY";
    }
    return "IMPOSSIBLE";
}

hid_t HDF::openFileHelper(const std::string& filename, HDF::Mode mode)  {
    hid_t file = -1;
    assert_msg(!filename.empty(), "Filename is empty.");
    std::ifstream test_access(filename);
    switch (mode) {
        case APPEND:
            if (test_access.good()) {  // exists
                file = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
            } else {
                file = H5Fcreate(filename.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
            }
            break;
        case DESTROY:
            file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
            break;
        case READONLY:
            assert_msg(test_access.good(), "To use readonly access, file '%s' must exist "
                                           "and be accessible.", filename);
            file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
            break;
    }
    assert_msg(file >= 0, "Opening file %s with mode %s failed with status %d. Are you sure "
                          "all potential intermediate directories exist?",
               filename, str(mode), file);
    return file;
}

void HDF::openGroup(std::string group_name) {
    assert_msg(H5Iis_valid(file), "File id %d invalid. Did you open a file before creating "
                                  "a group?", file);
    closeGroup();
    assert_msg(!group_name.empty(), "Group name must not be empty.");
    assert_msg(group_name.size() == 1 || group_name.back() != '/',
               "Group name must not end with a '/' unless it's root, got '%s'.", group_name);

    if (group_name[0] == '/') {
        group_name_ = "";
    } else {
        group_name = "/" + group_name;
        if (group_name_ == "/") group_name_ = "";
    }
    std::string::size_type idx = 0;
    do {
        auto new_idx = group_name.find('/', idx + 1);
        assert_msg(new_idx - idx > 1, "Two consecutive '/' not allowed in group name '%s'.",
                   group_name);
        if (new_idx == std::string::npos) { new_idx = group_name.size(); }
        group_name_ += group_name.substr(idx, new_idx - idx);
        idx = new_idx;

        // Root group treated separately to ensure compatibility with 1.8 version.
        // See: https://support.hdfgroup.org/HDF5/doc/RM/RM_H5L.html#Link-Exists
        if (group_name_ != "/" && !H5Lexists(file, group_name_.c_str(), H5P_DEFAULT)) {
            group = H5Gcreate(file, group_name_.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            assert_msg(group >= 0, "Failed creating group '%s' in file '%s' with status '%d'.",
                       group_name_, filename_, group);
            H5Gclose(group);
        }
    } while (idx != group_name.size());
    group = H5Gopen(file, group_name_.c_str(), H5P_DEFAULT);
    assert_msg(group >= 0, "Failed opening group '%s' in file '%s' with status '%d'.",
               group_name_, filename_, group);
}

/// Specialization for `std::string`. @sa HDF::readAttribute.
template <>
std::string HDF::readAttribute<std::string>(const std::string& attr_name) const {
    assert_msg(H5Iis_valid(group), "Group id %d invalid. Did you open a group before reading?",
               group);
    hid_t attr = H5Aopen(group, attr_name.c_str(), H5P_DEFAULT);
    assert_msg(attr >= 0, "Attribute '%s' could not be accessed in group '%s' in file '%s'.",
               attr_name, group_name_, filename_);

    hid_t type = H5Aget_type(attr);
    assert_msg(type >= 0, "Failed getting type of attribute '%s' in group '%s' in file '%s'.",
               attr_name, group_name_, filename_);

    H5T_str_t pad = H5Tget_strpad(type);
    assert_msg(pad >= 0, "Attribute '%s' in group '%s' in file '%s' is not a string attribute.",
               attr_name, group_name_, filename_);
    std::size_t size = H5Tget_size(type);
    std::vector<char> str;
    str.reserve(size);
    herr_t status = H5Aread(attr, type, str.data());
    assert_msg(status >= 0, "Failed reading attribute '%s' from group '%s' in file '%s'.",
               attr_name, groupName(), filename_);

    H5Tclose(type);
    H5Aclose(attr);
    return std::string(str.data(), size);
}

std::string HDF::readNullTerminatedStringAttribute(const std::string& attr_name) const {
    std::string s = readStringAttribute(attr_name);
    return std::string(s.c_str());
}

int HDF::readIntAttribute(const std::string& attr_name) const {
    return readAttribute<int>(attr_name);
}

bool HDF::readBoolAttribute(const std::string& attr_name) const {
    return readAttribute<bool>(attr_name);
}

double HDF::readDoubleAttribute(const std::string& attr_name) const {
    return readAttribute<double>(attr_name);
}

float HDF::readFloatAttribute(const std::string& attr_name) const {
    return readAttribute<float>(attr_name);
}

std::string HDF::readStringAttribute(const std::string& attr_name) const {
    return readAttribute<std::string>(attr_name);
}

void HDF::writeIntAttribute(const std::string& attr_name, int value, bool overwrite) const {
    writeAttribute(attr_name, value, H5T_NATIVE_INT, overwrite);
}

void HDF::writeBoolAttribute(const std::string& attr_name, bool value, bool overwrite) const {
    writeAttribute(attr_name, static_cast<hbool_t>(value), H5T_NATIVE_HBOOL, overwrite);
}

void HDF::writeDoubleAttribute(const std::string& attr_name, double value, bool overwrite) const {
    writeAttribute(attr_name, value, H5T_NATIVE_DOUBLE, overwrite);
}

void HDF::writeFloatAttribute(const std::string& attr_name, float value, bool overwrite) const {
    writeAttribute(attr_name, value, H5T_NATIVE_FLOAT, overwrite);
}

void HDF::writeStringAttribute(const std::string& attr_name, const std::string& value,
                               bool overwrite) const {
    assert_msg(H5Iis_valid(group), "Group id %d invalid. Did you open a group before writing?",
               group);
    assert_msg(!value.empty(), "Storing empty strings is not supported in HDF5.");
    if (H5Aexists(group, attr_name.c_str())) {
        if (!overwrite) {
            assert_msg(false, "Attribute '%s' in group '%s' in file '%s' already exists. To "
                              "overwrite its contents use parameter overwrite=true.",
                       attr_name, group_name_, filename_);
        }
        herr_t status = H5Adelete(group, attr_name.c_str());
        assert_msg(status >= 0, "Failed deleting existing attribute '%s' in group '%s' in "
                                "file '%s' before writing a new one.",
                   attr_name, group_name_, filename_);
    }

    hid_t space = H5Screate(H5S_SCALAR);
    hid_t type = H5Tcopy(H5T_C_S1);
    herr_t err = H5Tset_size(type, value.size());
    assert_msg(err >= 0, "Set size failed.");
    hid_t attr = H5Acreate2(group, attr_name.c_str(), type, space, H5P_DEFAULT, H5P_DEFAULT);
    assert_msg(attr >= 0, "Failed creating attribute '%s' in group '%s' in file '%s'.",
               attr_name, group_name_, filename_);

    herr_t status = H5Awrite(attr, type, value.data());
    assert_msg(status >= 0, "Failed writing attribute '%s' to group '%s' in file '%s'.",
               attr_name, group_name_, filename_);
    err = H5Tclose(type);
    assert_msg(err >= 0, "Type closing failed.");
    err = H5Sclose(space);
    assert_msg(err >= 0, "Space closing failed.");
    err = H5Aclose(attr);
    assert_msg(err >= 0, "Attribute closing failed.");
}

std::vector<int> HDF::readIntArray(const std::string& dataset_name) const {
    return read1DArray<int>(dataset_name);
}

std::vector<double> HDF::readDoubleArray(const std::string& dataset_name) const {
    return read1DArray<double>(dataset_name);
}

std::vector<float> HDF::readFloatArray(const std::string& dataset_name) const {
    return read1DArray<float>(dataset_name);
}

std::vector<std::vector<int>> HDF::readInt2DArray(const std::string& dataset_name) const {
    return read2DArray<int>(dataset_name);
}

std::vector<std::vector<double>> HDF::readDouble2DArray(const std::string& dataset_name) const {
    return read2DArray<double>(dataset_name);
}

std::vector<std::vector<float>> HDF::readFloat2DArray(const std::string& dataset_name) const {
    return read2DArray<float>(dataset_name);
}

std::vector<std::vector<std::vector<int>>>
HDF::readInt3DArray(const std::string& dataset_name) const {
    return read3DArray<int>(dataset_name);
}

std::vector<std::vector<std::vector<double>>>
HDF::readDouble3DArray(const std::string& dataset_name) const {
    return read3DArray<double>(dataset_name);
}

std::vector<std::vector<std::vector<float>>>
HDF::readFloat3DArray(const std::string& dataset_name) const {
    return read3DArray<float>(dataset_name);
}

void HDF::flush() const {
    if (!H5Iis_valid(file)) return;
    H5Fflush(file, H5F_SCOPE_GLOBAL);
    herr_t status = H5Fclose(file);
    assert_msg(status >= 0, "Flushing file '%s' failed with status %d.", filename_, status);
}

void HDF::closeGroup() const {
    if (!H5Iis_valid(group)) return;
    herr_t status = H5Gclose(group);
    assert_msg(status >= 0, "Closing group '%s' in file '%s' failed with status %d.",
               group_name_, filename_, status);
}

void HDF::closeFile() const {
    closeGroup();
    if (!H5Iis_valid(file)) return;
    herr_t status = H5Fclose(file);
    assert_msg(status >= 0, "Closing file '%s' failed with status %d.", filename_, status);
}

HDF HDF::atomic() const {
    assert_msg(!filename_.empty(), "Filename should not be empty.");
    assert_msg(!group_name_.empty(), "Group name should not be empty.");
    assert_msg(!H5Iis_valid(group), "Group '%s' already opened, using atomic makes no sense.",
               group_name_);
    assert_msg(!H5Iis_valid(file), "File '%s' already opened, using atomic makes no sense.",
               filename_);
    HDF hdf(filename_, APPEND);
    hdf.openGroup(group_name_);
    return hdf;  // file is closed in destructor
}

void HDF::setGroupName(const std::string& group_name) {
    closeGroup();
    assert_msg(!group_name.empty(), "Group name must not be empty");
    assert_msg(group_name.find("//") == std::string::npos, "Group name '%s' should not contain "
                                                           "two consecutive '/'.", group_name);
    if (group_name[0] == '/') {
        group_name_ = group_name;
    } else {
        group_name_ = '/' + group_name;
    }
    if (group_name.length() >= 2 && group_name.back() == '/') {
        group_name_.pop_back();
    }
}

herr_t HDF::memberIterateCallback(hid_t loc_id, const char* name, const H5L_info_t*, void* data) {
    auto* results = static_cast<Members*>(data);
    H5O_info_t infobuf;
#if H5_VERSION_GE(1, 12, 0)
    herr_t status = H5Oget_info_by_name(loc_id, name, &infobuf, H5O_INFO_ALL, H5P_DEFAULT);
#else
    herr_t status = H5Oget_info_by_name(loc_id, name, &infobuf, H5P_DEFAULT);
#endif
    assert_msg(status >= 0, "Could not get type of object '%s'.", name);
    switch (infobuf.type) {
        case H5O_TYPE_GROUP:
            results->groups.emplace_back(name);
            break;
        case H5O_TYPE_DATASET:
            results->datasets.emplace_back(name);
            break;
        case H5O_TYPE_NAMED_DATATYPE:
            results->datatypes.emplace_back(name);
            break;
        default:
            results->unknowns.emplace_back(name);
    }
    return 0;
}

HDF::Members HDF::members() const {
    assert_msg(H5Iis_valid(group), "Group id %d invalid. Did you open a group before "
                                   "iterating over it?", group);
    Members result;
    herr_t status = H5Literate(group, H5_INDEX_NAME, H5_ITER_INC, nullptr,
                               memberIterateCallback, static_cast<void*>(&result));
    assert_msg(status >= 0, "Error during iteration through group '%s'.", group_name_);
    return result;
}

}  // namespace mm
