#include <hdf5.h>

int main() {
    H5Fopen("tmp.h5", H5F_ACC_RDWR, H5P_DEFAULT);
    return 0;
}
