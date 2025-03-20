#ifndef AUX_MEDUSA_HEADER
#define AUX_MEDUSA_HEADER

//#include "../../medusa/include/medusa/Medusa.hpp"
#include "spmatrix_medusa.h"
#include <medusa/Medusa_fwd.hpp>
#include <Eigen/SparseCore>


// reduce size of stencils in domain d with idxs from size n_large to n_small
template<typename Domain>
void reduceStencil(Domain &d, int n_small, int n_large, mm::Range<int> idxs)
{
    mm::Range<int> rem_idxs = mm::Range<int>::seq(n_small,n_large);
    #pragma omp parallel for
    for(int i : idxs) {
        d.support(i).remove(rem_idxs); // only keep first n_small support nodes
    }
}

// convert matrix to block oseen format in medusa, gradient [3,4,5] and divergence [0,1,2] blocks
pspmatrixm matMM2MOseen(Eigen::SparseMatrix<double, Eigen::RowMajor> &mat, int N_ui, int N_p);

// adjust rhs according to reduceMat2
void adjustRHS2(mm::Range<int> idxs_ui, Eigen::VectorXd y, Eigen::VectorXd &rhs, Eigen::SparseMatrix<double, Eigen::RowMajor> Z, int N_u, int N_ub, int N_p);

// adjust rhs according to reduceMat3
void adjustRHS3(mm::Range<int>& idxs_ui, Eigen::VectorXd& y, Eigen::VectorXd &rhs, Eigen::SparseMatrix<double, Eigen::RowMajor>& Z, int N_u, int N_ub, int N_p);

// copy matrix B into block(rowStartA:B.rows()+rowStartA, colStartA:b.cols()+colStartA) of A
void insertBlock(Eigen::SparseMatrix<double, Eigen::RowMajor> &A, Eigen::SparseMatrix<double, Eigen::RowMajor>& B, int rowStartA, int colStartA);

// remove parts of matrix corresponding to boundary conditions 2D
void reduceMat2(Eigen::SparseMatrix<double, Eigen::RowMajor> &M, Eigen::SparseMatrix<double, Eigen::RowMajor> Z,
                Eigen::SparseMatrix<double, Eigen::RowMajor> B, int N_u, int N_p, int N_ui, int N_ub, Eigen::VectorXd constraint);

// remove parts of matrix corresponding to boundary conditions 3D
void reduceMat3(Eigen::SparseMatrix<double, Eigen::RowMajor> &M, Eigen::SparseMatrix<double, Eigen::RowMajor>& Z,
                Eigen::SparseMatrix<double, Eigen::RowMajor>& B, int N_u, int N_p, int N_ui, int N_ub, Eigen::VectorXd& constraint);

#endif
