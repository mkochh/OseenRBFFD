#include "aux_medusa.hpp"

pspmatrixm matMM2MOseen(Eigen::SparseMatrix<double, Eigen::RowMajor> &mat, int N_ui, int N_p)
{
    Eigen::SparseMatrix<double, Eigen::RowMajor> *B = new Eigen::SparseMatrix<double, Eigen::RowMajor>[6];
    // set gradient [3,4,5] and divergence [0,1,2] blocks
    B[0] = (Eigen::SparseMatrix<double, Eigen::RowMajor>)mat.block(3*N_ui,0,N_p+1,N_ui);
    B[3] = (Eigen::SparseMatrix<double, Eigen::RowMajor>)mat.block(0,3*N_ui,N_ui,N_p+1).transpose();
    B[1] = (Eigen::SparseMatrix<double, Eigen::RowMajor>)mat.block(3*N_ui,N_ui,N_p+1,N_ui);
    B[4] = (Eigen::SparseMatrix<double, Eigen::RowMajor>)mat.block(N_ui,3*N_ui,N_ui,N_p+1).transpose();
    B[2] = (Eigen::SparseMatrix<double, Eigen::RowMajor>)mat.block(3*N_ui,2*N_ui,N_p+1,N_ui);
    B[5] = (Eigen::SparseMatrix<double, Eigen::RowMajor>)mat.block(2*N_ui,3*N_ui,N_ui,N_p+1).transpose();

    Eigen::SparseMatrix<double, Eigen::RowMajor> *C = new Eigen::SparseMatrix<double, Eigen::RowMajor>;
    C[0] = (Eigen::SparseMatrix<double, Eigen::RowMajor>)mat.block(3*N_ui,3*N_ui,N_p+1,N_p+1);

    Eigen::SparseMatrix<double, Eigen::RowMajor> *A = new Eigen::SparseMatrix<double, Eigen::RowMajor>;
    A[0] = (Eigen::SparseMatrix<double, Eigen::RowMajor>)mat.block(0,0,N_ui,N_ui);
    
    pspmatrixm mat_oseen_m = new_spmatrixm(A, B, C); // stiffness matrix

    return mat_oseen_m;
}

void adjustRHS2(mm::Range<int> idxs_ui, Eigen::VectorXd y, Eigen::VectorXd &rhs, Eigen::SparseMatrix<double, Eigen::RowMajor> Z, int N_u, int N_ub, int N_p)
{
    mm::Range<int> idxs_ui_Nu = idxs_ui;
    mm::Range<int> idxs_ui_2Nu = idxs_ui;
    for(int i = 0; i < idxs_ui_Nu.size(); i++) {
        idxs_ui_Nu[i] += N_u;
        idxs_ui_2Nu[i] += 2*N_u;
    }
    Eigen::VectorXd u1b = y.segment(0, N_ub);
    Eigen::VectorXd u2b = y.segment(N_u, N_ub);
    int j = 0;
    for(int i : idxs_ui) {
        rhs(j) = y(i) - Z.row(i).segment(0,N_ub).dot(u1b);
        j++;
    }
    for(int i : idxs_ui_Nu) {
        rhs(j) = y(i) - Z.row(i).segment(N_u, N_ub).dot(u2b);
        j++;
    }
    for(int i : mm::Range<int>::seq(2*N_u,2*N_u+N_p)) {
        rhs(j) = y(i) - Z.row(i).segment(0, N_ub).dot(u1b) - Z.row(i).segment(N_u, N_ub).dot(u2b);
        j++;
    }
}

void adjustRHS3(mm::Range<int>& idxs_ui, Eigen::VectorXd& rhs_with_dirbnd, Eigen::VectorXd &rhs_without_dirbnd, Eigen::SparseMatrix<double, Eigen::RowMajor>& Z, int N_u, int N_ub, int N_p)
{
    mm::Range<int> idxs_ui_Nu = idxs_ui;
    mm::Range<int> idxs_ui_2Nu = idxs_ui;
    for(int i = 0; i < idxs_ui_Nu.size(); i++) {
        idxs_ui_Nu[i] += N_u;
        idxs_ui_2Nu[i] += 2*N_u;
    }
    Eigen::VectorXd u1b = rhs_with_dirbnd.segment(0, N_ub);
    Eigen::VectorXd u2b = rhs_with_dirbnd.segment(N_u, N_ub);
    Eigen::VectorXd u3b = rhs_with_dirbnd.segment(2*N_u, N_ub);
    int j = 0;
    for(int i : idxs_ui) {
        rhs_without_dirbnd(j) = rhs_with_dirbnd(i) - Z.row(i).segment(0,N_ub).dot(u1b);
        j++;
    }
    for(int i : idxs_ui_Nu) {
        rhs_without_dirbnd(j) = rhs_with_dirbnd(i) - Z.row(i).segment(N_u, N_ub).dot(u2b);
        j++;
    }
    for(int i : idxs_ui_2Nu) {
        rhs_without_dirbnd(j) = rhs_with_dirbnd(i) - Z.row(i).segment(2*N_u, N_ub).dot(u3b);
        j++;
    }
    for(int i : mm::Range<int>::seq(3*N_u,3*N_u+N_p)) {
        rhs_without_dirbnd(j) = rhs_with_dirbnd(i) - Z.row(i).segment(0, N_ub).dot(u1b) - Z.row(i).segment(N_u, N_ub).dot(u2b) - Z.row(i).segment(2*N_u, N_ub).dot(u3b);
        j++;
    }
}

void insertBlock(Eigen::SparseMatrix<double, Eigen::RowMajor> &A, Eigen::SparseMatrix<double, Eigen::RowMajor>& B, int rowStartA, int colStartA)
{
    #pragma omp parallel for
    for (int k=0; k<B.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(B,k); it; ++it) {
            int r = it.row();   // row index
            int c = it.col();   // col index
            A.insert(r+rowStartA,c+colStartA) = it.value();
        }
    }
}

void reduceMat2(Eigen::SparseMatrix<double, Eigen::RowMajor> &M, Eigen::SparseMatrix<double, Eigen::RowMajor> Z,
                Eigen::SparseMatrix<double, Eigen::RowMajor> B, int N_u, int N_p, int N_ui, int N_ub, Eigen::VectorXd constraint)
{
    B = Z.block(N_ub,N_ub,N_ui,N_ui);
    insertBlock(M,B,0,0);

    B = Z.block(N_u+N_ub,N_u+N_ub,N_ui,N_ui);
    insertBlock(M,B,N_ui,N_ui);

    B = Z.block(N_ub,2*N_u,N_ui,N_p);
    insertBlock(M,B,0,2*N_ui);

    B = Z.block(N_u+N_ub,2*N_u,N_ui,N_p);
    insertBlock(M,B,N_ui,2*N_ui);

    B = Z.block(2*N_u,N_ub,N_p,N_ui);
    insertBlock(M,B,2*N_ui,0);

    B = Z.block(2*N_u,N_u+N_ub,N_p,N_ui);
    insertBlock(M,B,2*N_ui,N_ui);

    for(int i = 2*N_ui; i < 2*N_ui+N_p; i++) {
        M.insert(2*N_ui+N_p,i) = constraint(i-2*N_ui);
        M.insert(i,2*N_ui+N_p) = constraint(i-2*N_ui);
    }
}

void reduceMat3(Eigen::SparseMatrix<double, Eigen::RowMajor> &M, Eigen::SparseMatrix<double, Eigen::RowMajor>& Z,
                Eigen::SparseMatrix<double, Eigen::RowMajor>& B, int N_u, int N_p, int N_ui, int N_ub, Eigen::VectorXd& constraint)
{
    B = Z.block(N_ub,N_ub,N_ui,N_ui);
    insertBlock(M,B,0,0);
    B = Z.block(N_u+N_ub,N_u+N_ub,N_ui,N_ui);
    insertBlock(M,B,N_ui,N_ui);
    B = Z.block(2*N_u+N_ub,2*N_u+N_ub,N_ui,N_ui);
    insertBlock(M,B,2*N_ui,2*N_ui);
    B = Z.block(N_ub,3*N_u,N_ui,N_p);
    insertBlock(M,B,0,3*N_ui);
    B = Z.block(N_u+N_ub,3*N_u,N_ui,N_p);
    insertBlock(M,B,N_ui,3*N_ui);
    B = Z.block(2*N_u+N_ub,3*N_u,N_ui,N_p);
    insertBlock(M,B,2*N_ui,3*N_ui);
    B = Z.block(3*N_u,N_ub,N_p,N_ui);
    insertBlock(M,B,3*N_ui,0);
    B = Z.block(3*N_u,N_u+N_ub,N_p,N_ui);
    insertBlock(M,B,3*N_ui,N_ui);
    B = Z.block(3*N_u,2*N_u+N_ub,N_p,N_ui);
    insertBlock(M,B,3*N_ui,2*N_ui);
    for(int i = 3*N_ui; i < 3*N_ui+N_p; i++) {
        M.insert(3*N_ui+N_p,i) = constraint(i-3*N_ui);
        M.insert(i,3*N_ui+N_p) = constraint(i-3*N_ui);
    }
}
