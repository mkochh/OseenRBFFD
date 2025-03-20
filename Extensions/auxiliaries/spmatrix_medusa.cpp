#include <stdio.h>

#include "spmatrix_medusa.h"

pspmatrixm new_spmatrixm(Eigen::SparseMatrix<double, Eigen::RowMajor> *A, Eigen::SparseMatrix<double, Eigen::RowMajor> *B, Eigen::SparseMatrix<double, Eigen::RowMajor> *C)
{

  pspmatrixm sp;

  sp = (pspmatrixm)malloc(sizeof(spmatrixm));

  sp->A = A;
  sp->B = B;
  sp->C = C;

  return sp;
}

void delete_spmatrixm(pspmatrixm sp)
{
  delete[] sp->B;
  delete sp->A;
  delete sp->C;
  free(sp);
}

Eigen::VectorXd addeval_spmatrixm_vector(double alpha, pspmatrixm sp, Eigen::VectorXd x, Eigen::VectorXd y)
{
  Eigen::VectorXd y1[3], y2;
  Eigen::VectorXd x1[3], x2;
  uint rows_A, rows_B, cols_A;
  uint d;

  rows_A = sp->A[0].rows();
  rows_B = sp->B[0].rows();
  cols_A = sp->A[0].cols();

  assert(x.size() == 3 * cols_A + rows_B);
  assert(y.size() == 3 * rows_A + rows_B);

  for (d = 0; d < 3; d++)
    x1[d] = x.segment(d * cols_A, cols_A);
  x2 = x.segment(3 * cols_A, rows_B);
  for (d = 0; d < 3; d++)
    y1[d] = y.segment(d * rows_A, rows_A);
  y2 = y.segment(3 * rows_A, rows_B);

  for (d = 0; d < 3; d++)
    y1[d] += alpha*sp->A[0]*x1[d];
  for (d = 0; d < 3; d++)
  {
    y1[d] += alpha*sp->B[2*d+1].transpose()*x2;
    y2 += alpha*sp->B[2*d]*x1[d];
  }
  y2 += alpha*sp->C[0]*x2;

  for (d = 0; d < 3; d++)
    y.segment(d* rows_A, rows_A) = y1[d];
  y.tail(rows_B) = y2;  

  return y;
}

/* void addevaltrans_spmatrixm_vector(double alpha, pspmatrixm sp, Eigen::VectorXd x, Eigen::VectorXd y)
{
  pavector y1[3], y2;
  pavector x1[3], x2;
  uint rows_A, rows_B, cols_A;
  uint d;

  rows_A = sp->A->rows;
  rows_B = sp->B[0]->rows;
  cols_A = sp->A->cols;

  for (d = 0; d < 3; d++)
    x1[d] = new_sub_avector((pavector)x, cols_A, d * cols_A);
  x2 = new_sub_avector((pavector)x, rows_B, 3 * cols_A);
  for (d = 0; d < 3; d++)
    y1[d] = new_sub_avector(y, rows_A, d * rows_A);
  y2 = new_sub_avector(y, rows_B, 3 * rows_A);

  for (d = 0; d < 3; d++)
    mvm_sparsematrix_avector(alpha, true, sp->A, x1[d], y1[d]);
  for (d = 0; d < 3; d++)
  {
    mvm_sparsematrix_avector(alpha, true, sp->B[d], x2, y1[d]);
    mvm_sparsematrix_avector(alpha, false, sp->B[d], x1[d], y2);
  }
  mvm_sparsematrix_avector(-1.0 * alpha, true, sp->C, x2, y2);
} */

/* void mvm_spmatrixm_vector(double alpha, bool trans, pspmatrixm sp, Eigen::VectorXd x, Eigen::VectorXd y)
{

  if (trans)
    addevaltrans_spmatrixm_vector(alpha, sp, x, y);
  else
    addeval_spmatrixm_vector(alpha, sp, x, y);
} */