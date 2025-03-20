#ifndef SET_WEIGHTS_HEADER
#define SET_WEIGHTS_HEADER

//#include "../../medusa/include/medusa/Medusa.hpp"
#include "rbffd_discretization.hpp"
#include <medusa/Medusa_fwd.hpp>
#include "../../medusa/include/medusa/bits/approximations/RBFFD.hpp"
#include "../../medusa/include/medusa/bits/approximations/Polyharmonic.hpp"
#include <Eigen/SparseCore>

namespace mm {

template<typename Operator1, typename Operator2>
void setWeightsCD3(Operator1& op_conv, Operator2& op_lap, const OseenDiscretizationBetter& dc)
{
    #pragma omp parallel for
    for (int i : dc.idxs_ui) {
        -dc.nu * op_lap.lap(i) = dc.rhs_vec3d[i];
        op_conv.grad(i,dc.conv_vec3d[i]) = 0;
    }
    #pragma omp parallel for
    for (int i : dc.d_u.boundary()) {
        op_lap.value(i) = dc.dir_bnd_vec3d[i];
    }
}

template<typename Operator>
void setWeightsDIV3(Operator& op_div, int N_u, Range<int>& idxs_p)
{
    op_div.setRowOffset(2*N_u);
    #pragma omp parallel for
    for (int i : idxs_p) {
        op_div.der1(i,0) = 0;
    }
    op_div.setColOffset(N_u);
    #pragma omp parallel for
    for (int i : idxs_p) {
        op_div.der1(i,1) = 0;
    }
    op_div.setColOffset(2*N_u);
    #pragma omp parallel for
    for (int i : idxs_p) {
        op_div.der1(i,2) = 0;
    }
}

template<typename Operator>
void setWeightsGRAD3(Operator& op_grad, int N_u, Range<int>& idxs_ui)
{
    op_grad.setColOffset(2*N_u);
    #pragma omp parallel for
    for (int i : idxs_ui) {
        op_grad.der1(i,0) = 0;
    }
    op_grad.setRowOffset(N_u);
    #pragma omp parallel for
    for (int i : idxs_ui) {
        op_grad.der1(i,1) = 0;
    }
    op_grad.setRowOffset(2*N_u);
    #pragma omp parallel for
    for (int i : idxs_ui) {
        op_grad.der1(i,2) = 0;
    }
}

void setPressureConstraint3(Eigen::SparseMatrix<double, Eigen::RowMajor> &mat_with_dirbnd, int N_u, int N_p, Eigen::VectorXd &constraint)
{
    for(int i = 3*N_u; i < 3*N_u+N_p; i++) {
        mat_with_dirbnd.insert(3*N_u+N_p,i) = constraint(i-3*N_u);
        mat_with_dirbnd.insert(i,3*N_u+N_p) = constraint(i-3*N_u);
    }
}

Eigen::VectorXd setWeightsOseen_versatile_better(Eigen::SparseMatrix<double, Eigen::RowMajor> &mat_with_dirbnd, OseenDiscretizationBetter &dc)
{
    int l_lap = dc.l[dc.poly_lap], l_grad = dc.l[dc.poly_grad], l_conv = dc.l[dc.poly_conv], l_div = dc.l[dc.poly_div];
    int k_lap = dc.k[dc.poly_lap], k_grad = dc.k[dc.poly_grad], k_conv = dc.k[dc.poly_conv], k_div = dc.k[dc.poly_div];
    int n1 = dc.n[dc.poly_grad], n2 = std::max(dc.n[dc.poly_lap],dc.n[dc.poly_conv]), n3 = dc.n[dc.poly_div];
    Eigen::VectorXi reserve_vector(3*dc.N_u+dc.N_p+1);
    reserve_vector << Eigen::VectorXi::Constant(dc.N_ub,1), Eigen::VectorXi::Constant(dc.N_ui,n2+n1), 
                        Eigen::VectorXi::Constant(dc.N_ub,1), Eigen::VectorXi::Constant(dc.N_ui,n2+n1), 
                        Eigen::VectorXi::Constant(dc.N_ub,1), Eigen::VectorXi::Constant(dc.N_ui,n2+n1), 
                        Eigen::VectorXi::Constant(dc.N_p,3*n3+1), Eigen::VectorXi::Constant(1,dc.N_p);
    mat_with_dirbnd.reserve(reserve_vector);
    Eigen::VectorXd rhs_with_dirbnd(3*dc.N_u+dc.N_p+1); rhs_with_dirbnd.setZero();
    Eigen::VectorXd rhs_dummy(3*dc.N_u+dc.N_p+1); rhs_dummy.setZero(); // dummy vector, discarded later
    
    Polyharmonic<double, -1> ph_lap(k_lap);
    RBFFD<Polyharmonic<double, -1>, Vec3d, ScaleToClosest, Eigen::PartialPivLU<Eigen::MatrixXd>> approx_lap(ph_lap, l_lap);
    auto storage_lap = dc.d_u.computeShapes<sh::lap>(approx_lap, dc.idxs_ui);
    auto op_lap = storage_lap.implicitVectorOperators(mat_with_dirbnd, rhs_with_dirbnd);

    Polyharmonic<double, -1> ph_conv(k_conv);
    RBFFD<Polyharmonic<double, -1>, Vec3d, ScaleToClosest, Eigen::PartialPivLU<Eigen::MatrixXd>> approx_conv(ph_conv, l_conv);

    auto storage_conv = dc.d_conv.computeShapes<sh::d1>(approx_conv, dc.idxs_ui);
    auto op_conv = storage_conv.implicitVectorOperators(mat_with_dirbnd,rhs_with_dirbnd);
    
    Polyharmonic<double, -1> ph_grad(k_grad);
    RBFFD<Polyharmonic<double, -1>, Vec3d, ScaleToFarthest, Eigen::PartialPivLU<Eigen::MatrixXd>> approx_grad(ph_grad, l_grad);
    
    auto storage_grad = dc.d_grad.computeShapes<sh::d1>(approx_grad, dc.idxs_ui);
    auto op_grad = storage_grad.implicitOperators(mat_with_dirbnd, rhs_dummy);

    Polyharmonic<double, -1> ph_div(k_div);
    RBFFD<Polyharmonic<double, -1>, Vec3d, ScaleToFarthest, Eigen::PartialPivLU<Eigen::MatrixXd>> approx_div(ph_div, l_div);

    auto storage_div = dc.d_div.computeShapes<sh::d1>(approx_div, dc.idxs_p);
    auto op_div = storage_div.implicitOperators(mat_with_dirbnd, rhs_dummy);

    setWeightsCD3(op_conv, op_lap, dc);
    setWeightsDIV3(op_div, dc.N_u, dc.idxs_p);
    setWeightsGRAD3(op_grad, dc.N_u, dc.idxs_ui);
    setPressureConstraint3(mat_with_dirbnd, dc.N_u, dc.N_p, dc.constraint);

    return rhs_with_dirbnd;
}

// remove part in matrix corresponding to boundary conditions
Eigen::VectorXd eliminateDirichlet(Eigen::SparseMatrix<double, Eigen::RowMajor> &mat_without_dirbnd, Eigen::SparseMatrix<double, Eigen::RowMajor> &mat_with_dirbnd, 
                                    OseenDiscretizationBetter &dc, Eigen::VectorXd &rhs_with_dirbnd)
{
    int n1 = dc.n[dc.poly_grad], n2 = std::max(dc.n[dc.poly_lap],dc.n[dc.poly_conv]), n3 = dc.n[dc.poly_div];
    Eigen::VectorXi reserve_vector_mat(3*dc.N_ui+dc.N_p+1);
    reserve_vector_mat << Eigen::VectorXi::Constant(dc.N_ui, n2+n1), 
                        Eigen::VectorXi::Constant(dc.N_ui, n2+n1), 
                        Eigen::VectorXi::Constant(dc.N_ui, n2+n1), 
                        Eigen::VectorXi::Constant(dc.N_p, 3*n3+1), Eigen::VectorXi::Constant(1, dc.N_p);
    mat_without_dirbnd.reserve(reserve_vector_mat);
    Eigen::SparseMatrix<double, Eigen::RowMajor> T; // temporarily stores blocks of mat_with_dirbnd in function reduceMat, not needed afterwards
    reduceMat3(mat_without_dirbnd, mat_with_dirbnd, T, dc.N_u, dc.N_p, dc.N_ui, dc.N_ub, dc.constraint);
    mat_without_dirbnd.makeCompressed();

    Eigen::VectorXd rhs_without_dirbnd(3*dc.N_ui+dc.N_p+1); rhs_without_dirbnd.setZero();
    adjustRHS3(dc.idxs_ui, rhs_with_dirbnd, rhs_without_dirbnd, mat_with_dirbnd, dc.N_u, dc.N_ub, dc.N_p); // adjust rhs accordingly

    return rhs_without_dirbnd;
}

} // namespace mm

#endif
