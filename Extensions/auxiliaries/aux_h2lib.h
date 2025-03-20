#ifndef AUX_H2LIB_HEADER
#define AUX_H2LIB_HEADER

#include "../../H2Lib/h2lib.h"
//#include "../../medusa/include/medusa/Medusa.hpp"
#include "../sparse/spmatrix.h"
#include <medusa/Medusa_fwd.hpp>
#include <Eigen/SparseCore>


uint* idxMM2H(mm::indexes_t idxs_mm);

pavector vecMM2H(Eigen::VectorXd &vec_eigen);

Eigen::VectorXd vecH2MM(pavector vec_h2);

psparsematrix matMM2H(Eigen::SparseMatrix<double, Eigen::RowMajor> &mat);

// convert matrix to block oseen format in h2lib, gradient [3,4,5] and divergence [0,1,2] blocks
pspmatrix matMM2HOseen(Eigen::SparseMatrix<double, Eigen::RowMajor> &mat, int N_ui, int N_p);

// construct clustergeometry for H-Matrix procedure based on medusa domain
pclustergeometry build_clustergeometry_medusa(mm::DomainDiscretization<mm::Vec3d> domain, mm::Range<int> idx);

#endif
