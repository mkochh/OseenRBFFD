#ifndef SPMATRIX_MEDUSA_HEADER
#define SPMATRIX_MEDUSA_HEADER

//#include "../../medusa/include/medusa/Medusa.hpp"
#include <medusa/Medusa_fwd.hpp>
#include <Eigen/SparseCore>

/* Representation of a saddle point matrix, i.e. of a 2x2 block matrix
   with (1,1) block A, (2,1) block B, (2,1) block B^T, (2,2) block 0 */
struct _spmatrix_medusa
{
  Eigen::SparseMatrix<double, Eigen::RowMajor> *A;  // (1,1) block of the matrix
  Eigen::SparseMatrix<double, Eigen::RowMajor> *B; // (2,1) block of the matrix, and transposed of the (1,2) block
  Eigen::SparseMatrix<double, Eigen::RowMajor> *C;  // Stabilization part
};

typedef struct _spmatrix_medusa spmatrixm;
typedef spmatrixm *pspmatrixm;
typedef const pspmatrixm pcspmatrixm;

/* Constructor and destructor for a saddale point matrix */
pspmatrixm new_spmatrixm(Eigen::SparseMatrix<double, Eigen::RowMajor> *A, Eigen::SparseMatrix<double, Eigen::RowMajor> *B, Eigen::SparseMatrix<double, Eigen::RowMajor> *C);

void delete_spmatrixm(pspmatrixm sp);

/* Matrix vector multiplication for saddle point matrices */
Eigen::VectorXd addeval_spmatrixm_vector(double alpha, pspmatrixm sp, Eigen::VectorXd x, Eigen::VectorXd y);

#endif