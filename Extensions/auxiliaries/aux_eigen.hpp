#ifndef AUX_EIGEN_HEADER
#define AUX_EIGEN_HEADER

#include "spmatrix_medusa.h"
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

// Eigen class BiCGstab doesn't allow system matrix and preconditioner matrix
// to be different, this is a workaround using the internal implementation of Eigen
template<typename MatrixType, typename Preconditioner>
bool bicgstab(const MatrixType& mat, const Eigen::VectorXd& rhs,
              Eigen::VectorXd& x, const Preconditioner& precond, int& iters, double& tol_error)
{
    double tol = tol_error;
    int maxIters = iters;

    int n = mat.cols();
    Eigen::VectorXd r = rhs - mat * x;
    Eigen::VectorXd r0 = r;

    double r0_sqnorm = r0.squaredNorm();
    double rhs_sqnorm = rhs.squaredNorm();
    if(rhs_sqnorm == 0) {
        x.setZero();
        iters = 0;
        tol_error = 0;
        return true;
    }

    double rho = 1, alpha = 1, w = 1;

    Eigen::VectorXd v = Eigen::VectorXd::Zero(n), p = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd y(n), z(n), kt(n), ks(n), s(n), t(n);

    double tol2 = tol * tol * rhs_sqnorm;
    double eps2 = Eigen::NumTraits<double>::epsilon() * Eigen::NumTraits<double>::epsilon();
    int i = 0;
    // int restarts = 0;

    while ( r.squaredNorm() > tol2 && i<maxIters )
    {
        double rho_old = rho;

        rho = r0.dot(r);
        if (std::abs(rho) < eps2*r0_sqnorm)
        {
            // The new residual vector became too orthogonal to the arbitrarily chosen direction r0
            // Let's restart with a new r0:
            r  = rhs - mat * x;
            r0 = r;
            rho = r0_sqnorm = r.squaredNorm();
            /* if(restarts++ == 0)
            {
                i = 0;
            } */
        }
        double beta = (rho/rho_old) * (alpha / w);
        p = r + beta * (p - w * v);

        y = precond.solve(p);

        v.noalias() = mat * y;

        alpha = rho / r0.dot(v);
        s = r - alpha * v;

        z = precond.solve(s);
        t.noalias() = mat * z;

        double tmp = t.squaredNorm();
        if(tmp>double(0)) {
            w = t.dot(s) / tmp;
        } else {
            w = double(0);
        }
        x += alpha * y + w * z;
        r = s - w * t;

        /* std::cout << "iter = " << i << "\tr.squaredNorm() = " << r.squaredNorm() 
        << "\nrho = " << rho 
        << "\nr0.dot(v) = " << r0.dot(v)
        << "\ntmp = " << tmp
        << "\n" << std::endl; */

        ++i;
    }
    tol_error = std::sqrt(r.squaredNorm()/rhs_sqnorm);
    iters = i;
    if(iters>=maxIters || tol_error > tol) {
        return false;
    } else {
        return true;
    }
}

// ILUTCustom extends Eigen::IncompleteLUT
// allows to get access to permutations m_P/m_Pinv and factorization m_lu
// allows for custom permutations
template <typename _Scalar, typename _StorageIndex = int>
class ILUTCustom : public Eigen::IncompleteLUT<_Scalar, _StorageIndex>
{
    public:
        enum ordering {
            // RcmForDirect = 0,
            // RcmForBlocked = 1,
            AMD = 2,
            // Metis = 3
        };

    public:
        typedef Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic,_StorageIndex> PermType;

        typename Eigen::IncompleteLUT<_Scalar, _StorageIndex>::FactorType getLU() { return this->m_lu; }
        PermType getPerm() { return this->m_P; }
        PermType getPermInv() { return this->m_Pinv; }

        // could overload compute/analyzePattern to use different orderings
        // or maybe adding ordering method as template parameter could work
        template<typename MatrixType>
        void analyzePattern(const MatrixType& amat, int ordering, int N_ui);

        template<typename MatrixType>
        ILUTCustom& compute(const MatrixType& amat, int ordering = AMD, int N_ui = 0)
        {
            analyzePattern(amat, ordering, N_ui); 
            this->factorize(amat);
            return *this;
        }
};

template <typename Scalar, typename StorageIndex>
template<typename _MatrixType>
void ILUTCustom<Scalar,StorageIndex>::analyzePattern(const _MatrixType& amat, int ordering, int N_ui)
{
    assert(amat.isCompressed() && "Matrix has to be in compressed mode for orderings to work");
    // for direct ILUT amat will be whole matrix including all velocity blocks and pressure, have to take only one velocity block
    // simply applying RCM to whole matrix leads to bad reordering
    // if (ordering == RcmForDirect) {
    //     Eigen::SparseMatrix<double, Eigen::RowMajor> mat_block;
    //     mat_block = amat.block(0, 0, N_ui, N_ui);
    //     assert(mat_block.isCompressed() && "Block of compressed matrix is not compressed");
    //     Eigen::SparseMatrix<double, Eigen::RowMajor> temp = mat_block.transpose();
    //     Eigen::SparseMatrix<double, Eigen::RowMajor> mat_block_sym = mat_block + temp;

    //     amgcl::adapter::reorder<> perm(mat_block_sym);

    //     Eigen::VectorXi vec(N_ui), perm_vec_block_temp(N_ui), perm_vec_block(N_ui), perm_vec_full(amat.rows());
    //     for (int i = 0; i < N_ui; i++)
    //         vec(i) = i;

    //     perm.forward(vec, perm_vec_block_temp); // get cuthill mckee
    //     for (int i = 0; i < N_ui; i++)
    //         perm_vec_block(i) = perm_vec_block_temp(N_ui-i-1); // reverse cuthill mckee

    //     for (int i = 0; i < N_ui; i++) { perm_vec_full(i) = perm_vec_block(i); }
    //     for (int i = 0; i < N_ui; i++) { perm_vec_full(i+N_ui) = perm_vec_block(i) + N_ui; }
    //     for (int i = 0; i < N_ui; i++) { perm_vec_full(i+2*N_ui) = perm_vec_block(i) + 2*N_ui; }
    //     for (int i = 3*N_ui; i < amat.rows(); i++) { perm_vec_full(i) = i; }

    //     this->m_P = perm_vec_full.asPermutation();
    // }

    // for blockwise prcd ILUT amat will be only velocity block or schur block can work with those directly
    // if (ordering == RcmForBlocked) {
    //     Eigen::SparseMatrix<double, Eigen::RowMajor> temp = amat.transpose();
    //     Eigen::SparseMatrix<double, Eigen::RowMajor> mat_sym = amat + temp;

    //     amgcl::adapter::reorder<> perm(mat_sym);

    //     Eigen::VectorXi vec(amat.rows()), perm_vec_temp(amat.rows()), perm_vec(amat.rows());
    //     for (int i = 0; i < amat.rows(); i++)
    //         vec(i) = i;

    //     perm.forward(vec, perm_vec_temp); // get cuthill mckee
    //     for (int i = 0; i < amat.rows(); i++)
    //         perm_vec(i) = perm_vec_temp(amat.rows()-i-1); // reverse cuthill mckee

    //     this->m_P = perm_vec.asPermutation();
    // }

    if (ordering == AMD) {
        Eigen::SparseMatrix<Scalar,Eigen::ColMajor, StorageIndex> mat1 = amat;
        Eigen::AMDOrdering<StorageIndex> ordering;
        ordering(mat1,this->m_P); // symmetrizing of the matrix happens internally
    }

    /* if (ordering == Metis) {
        Eigen::MetisOrdering<StorageIndex> ordering;
        ordering(amat,this->m_P); // symmetrizing of the matrix happens internally
    } */
  
    this->m_Pinv  = this->m_P.inverse(); // cache the inverse permutation
    this->m_analysisIsOk = true;
    this->m_factorizationIsOk = false;
    this->m_isInitialized = true;
}

// compute ILUT BlockTriangularPreconditioner for saddle point matrix
template <typename _Scalar, typename _StorageIndex = int>
class OseenILUT : public Eigen::SparseSolverBase< OseenILUT<_Scalar, _StorageIndex> >
{
    public:
        typedef _Scalar Scalar;
        typedef _StorageIndex StorageIndex;
        typedef typename Eigen::NumTraits<Scalar>::Real RealScalar;
        typedef Eigen::SparseMatrix<Scalar,Eigen::RowMajor,StorageIndex> FactorType;
        typedef Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic,StorageIndex> PermType;
    
        enum {
          ColsAtCompileTime = Eigen::Dynamic,
          MaxColsAtCompileTime = Eigen::Dynamic
        };

    public:
    
        OseenILUT()
        : m_droptolA(Eigen::NumTraits<Scalar>::dummy_precision()), m_fillfactorA(10),
         m_droptolS(Eigen::NumTraits<Scalar>::dummy_precision()), m_fillfactorS(10)
        {}

        OseenILUT& compute(const pspmatrixm K, int ordering = ILUTCustom<Scalar>::ordering::AMD)
        {
            m_K = K;
            precondA.setDroptol(m_droptolA);
            precondA.setFillfactor(m_fillfactorA);
            precondA.compute(K->A[0], ordering);
            FactorType lu = precondA.getLU();
            PermType P = precondA.getPerm();
            PermType Pinv = precondA.getPermInv();

            // compute schur complement S
            Eigen::SparseMatrix<double, Eigen::RowMajor> S = K->C[0];
            // S -= (U^-T*(B_{2i}*P)^T)^T * (L^-1*(P^-1*B_{2i+1}^T))
            for (int i = 0; i < 3; i++)
            {
                Eigen::SparseMatrix<double, Eigen::ColMajor> B0 = (K->B[i]*P).transpose();
                Eigen::SparseMatrix<double, Eigen::ColMajor> B1 = (Pinv*K->B[i+3].transpose());
                computeSchurcomplementDenseChunk(S, B0, B1, lu, 1000);
            }
            precondS.setDroptol(m_droptolS);
            precondS.setFillfactor(m_fillfactorS);
            precondS.compute(S);

            return *this;
        }

        void setDroptolA(const RealScalar& droptol) { this->m_droptolA = droptol; }
        void setDroptolS(const RealScalar& droptol) { this->m_droptolS = droptol; }
        void setFillfactorA(int fillfactor) { this->m_fillfactorA = fillfactor; }
        void setFillfactorS(int fillfactor) { this->m_fillfactorS = fillfactor; }

        int cols() const { return 3*precondA.cols() + precondS.rows(); }
        int rows() const { return 3*precondA.rows() + precondS.rows(); }

        // overload solve, see Eigen documentation IdentityPreconditioner and DiagonalPreconditioner
        template<typename Rhs> inline const Eigen::Solve<OseenILUT, Rhs>
        solve(const Eigen::MatrixBase<Rhs>& b) const {
            return Eigen::Solve<OseenILUT, Rhs>(*this, b.derived());
        }

        template<typename Rhs, typename Dest>
        void _solve_impl(const Rhs& b, Dest& x) const
        {
            Rhs b1[3], b2;
            int Nui = precondA.rows();
            int Np = precondS.rows();
            for (uint i = 0; i < 3; i++) {
                b1[i] = b.segment(i*Nui,Nui);
                x.segment(i*Nui,Nui) = precondA.solve(b1[i]);
            }
            b2 = b.tail(Np);
            for (uint i = 0; i < 3; i++) {
                b2 -= m_K->B[i]*x.segment(i*Nui,Nui);
            }
            x.tail(Np) = precondS.solve(b2);
        }

    protected:

        RealScalar m_droptolA;
        int m_fillfactorA;
        RealScalar m_droptolS;
        int m_fillfactorS;
        ILUTCustom<Scalar> precondA;
        ILUTCustom<Scalar> precondS;
        pspmatrixm m_K;

        // compute part of a dense Schur complement S, given an LU factorization and rectangular matrices B0 and B1
        // auxiliary function needed for class OseenILUT
        // split matrices B0, B1 into chunks of size th and compute/update schurcomplement for each chunk
        // this limits storage uses because dense matrices have to be used
        template<typename MatrixTypeS, typename MatrixTypeB0, typename MatrixTypeB1, typename MatrixTypeLU>
        void computeSchurcomplementDenseChunk(MatrixTypeS& S, MatrixTypeB0& B0, MatrixTypeB1& B1, const MatrixTypeLU& LU, int th) {
            for (int j = 0; j < B0.cols(); j+=th) {
                int n = (th < B0.cols()-j) ? th : B0.cols()-j; // max(th, B0.cols-j) <- chunk size
                Eigen::MatrixXd b0 = B0.middleCols(j,n);
                Eigen::MatrixXd b1 = B1.middleCols(j,n);
                B0.middleCols(j,n) = LU.transpose().template triangularView<Eigen::Lower>().solve(b0).sparseView();
                B1.middleCols(j,n) = LU.template triangularView<Eigen::UnitLower>().solve(b1).sparseView();
            }
            S -= (MatrixTypeS)(B0.transpose() * B1);
        }
};

#endif
