#include "aux_h2lib.h"

uint* idxMM2H(mm::indexes_t idxs_mm)
{
    int N = idxs_mm.size();
    uint* idxs_h2 = (uint*)malloc(N*sizeof(uint));
    #pragma omp parallel for
    for (int i = 0; i < N; i++)
    {
        idxs_h2[i] = idxs_mm[i];
    }

    return idxs_h2;
}

pavector vecMM2H(Eigen::VectorXd &vec_eigen)
{
    int N = vec_eigen.size();
    pavector vec_h2 = new_avector(N);

    #pragma omp parallel for
    for (int i = 0; i < N; i++)
    {
        vec_h2->v[i] = vec_eigen[i];
    }
    return vec_h2;
}

Eigen::VectorXd vecH2MM(pavector vec_h)
{
    int N = vec_h->dim;
    Eigen::VectorXd vec_eigen(N);

    #pragma omp parallel for
    for (int i = 0; i < N; i++)
    {
        vec_eigen[i] = vec_h->v[i];
    }
    return vec_eigen;
}

psparsematrix matMM2H(Eigen::SparseMatrix<double, Eigen::RowMajor> &mat)
{
    uint rows = mat.rows();
    uint cols = mat.cols();
    psparsepattern sp = new_sparsepattern(rows, cols);
    #pragma omp parallel for
    for (int k=0; k<mat.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(mat,k); it; ++it) {
            addnz_sparsepattern(sp, it.row(), it.col());
        }
    }

    psparsematrix mat_h2 = new_zero_sparsematrix(sp);
    #pragma omp parallel for
    for (int k=0; k<mat.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(mat,k); it; ++it) {
            setentry_sparsematrix(mat_h2, it.row(), it.col(), it.value());
        }
    }

    del_sparsepattern(sp);

    return mat_h2;
}

pspmatrix matMM2HOseen(Eigen::SparseMatrix<double, Eigen::RowMajor> &mat, int N_ui, int N_p)
{
    psparsematrix *B = (psparsematrix *)allocmem(6 * sizeof(psparsematrix));
    Eigen::SparseMatrix<double, Eigen::RowMajor> tempMat;
    // set gradient [3,4,5] and divergence [0,1,2] blocks
    tempMat = mat.block(3*N_ui,0,N_p+1,N_ui);
    B[0] = matMM2H(tempMat);
    tempMat = mat.block(0,3*N_ui,N_ui,N_p+1).transpose();
    B[3] = matMM2H(tempMat);
    tempMat = mat.block(3*N_ui,N_ui,N_p+1,N_ui);
    B[1] = matMM2H(tempMat);
    tempMat = mat.block(N_ui,3*N_ui,N_ui,N_p+1).transpose();
    B[4] = matMM2H(tempMat);
    tempMat = mat.block(3*N_ui,2*N_ui,N_p+1,N_ui);
    B[2] = matMM2H(tempMat);
    tempMat = mat.block(2*N_ui,3*N_ui,N_ui,N_p+1).transpose();
    B[5] = matMM2H(tempMat);

    // set set lower right block
    tempMat = mat.block(3*N_ui,3*N_ui,N_p+1,N_p+1);
    psparsematrix C = matMM2H(tempMat);

    // set velocity block
    tempMat = mat.block(0,0,N_ui,N_ui);
    psparsematrix A = matMM2H(tempMat);
    
    pspmatrix mat_oseen_h2 = new_spmatrix(A, B, C); // stiffness matrix

    return mat_oseen_h2;
}

pclustergeometry build_clustergeometry_medusa(mm::DomainDiscretization<mm::Vec3d> domain, mm::Range<int> idx)
{
    int N = idx.size();
    int dim = domain.dim;
    pclustergeometry cg = new_clustergeometry(dim, N);
    #pragma omp parallel for
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            cg->x[i][j] = domain.pos(idx[i],j);
            cg->smin[i][j] = domain.pos(idx[i],j);
            cg->smax[i][j] = domain.pos(idx[i],j);
        }
    }

    return cg;
}
