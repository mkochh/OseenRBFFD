#ifndef HLU_HEADER
#define HLU_HEADER

#include "../../H2Lib/h2lib.h"
//#include "../../medusa/include/medusa/Medusa.hpp"
#include <medusa/Medusa_fwd.hpp>
#include "prec.h"

/**
 * @brief Collection of parameters to build up the H-Matrices
 *        for an Block_HLU_Prec preconditioner
 *
 */
struct HLU_Options
{

  pspmatrix K; // Saddle point system matrix

  uint clfu; // Maximal number of points in a leaf cluster (velocity)
  uint clfp; // Maximal number of points in a leaf cluster (pressure)

  clustermode modeu; // Algorithm for clustering of the velocity index set
  clustermode modep; // Algorithm for clustering of the pressure index set

  // Admissibility parameters
  void *eta_vel;
  void *eta_grad;
  void *eta_schur;

  admissible admu;  // Admissibility condition for the Laplace/Convection part
  admissible admpu; // Admissibility condition for the grad/div matrix
  admissible admp;  // Admissibility condition for the Schurcomplement

  ptruncmode tm;  // Truncation mode to use
  real eps_vel;   // Truncation accuracy for the velocity block
  real eps_mul;   // Truncation accuracy for the Schur complement computation
  real eps_schur; // Truncation accuracy for the Schur complement

  /**
   * @brief Construct a new hlr options object
   *
   * @param K Saddle point system matrix
   * @param clfu Maximal number of points in a leaf cluster (velocity)
   * @param clfp Maximal number of points in a leaf cluster (pressure)
   * @param modeu Algorithm for clustering of the velocity index set
   * @param modep Algorithm for clustering of the pressure index set
   * @param eta_vel Admissibility parameter for the velocity matrix
   * @param eta_grad Admissibility paramter for the gradient matrix
   * @param eta_schur Admissibility paramter for the Schur complement
   * @param admu Admissibility condition for the velocity block
   * @param admpu Admissibility condition for the grad/div block
   * @param admp Admissibility condition for the Schurcomplement
   * @param tm Truncation mode to use
   * @param eps Truncation accuracy
   */

  HLU_Options(pspmatrix K, uint clfu, uint clfp,
              clustermode modeu, clustermode modep,
              void *eta_vel, void *eta_grad, void *eta_schur,
              admissible admu, admissible admpu, admissible admp, ptruncmode tm,
              real eps_vel, real eps_mul, real eps_schur);
  /**
   * @brief Destroy the hlr options object
   *
   */
  ~HLU_Options(){};
};

struct Block_HLU_Times
{
  // Timings for the cluster tree generation
  real pressure_cluster;
  real velocity_cluster;

  real velocity_lu;

  real grad_upper_solve;
  real grad_lower_solve;
  real grad_schur_multiplication;
  real schur_computation;

  real schur_lu;

  real velocity_compression;
  real gradient_compression;
};

enum BlockPrcdType
{
  BLOCK_DIAGONAL,
  BLOCK_TRIANGULAR
};

class Block_HLU_Prcd : public Preconditioner
{
public:
  phmatrix A;
  phmatrix S;
  psparsematrix *B;
  psparsematrix C;

  pcluster root_velocity, root_pressure;
  pblock block_velocity, block_pressure, block_grad;
  uint *idxu, *idxp;

  BlockPrcdType prcd_type;

  /**
   * @brief Construct a new Block_HLU_Prec preconditioner and measure
   * construction timings
   */

  Block_HLU_Prcd(mm::DomainDiscretization<mm::Vec3d> d_u,
                          mm::DomainDiscretization<mm::Vec3d> d_div,
                          mm::Range<int> idxp_supp,
                          HLU_Options *opt,
                          Block_HLU_Times *times,
                          psparsematrix symA,
                          BlockPrcdType prcd_type = BLOCK_DIAGONAL);

  ~Block_HLU_Prcd();

  /**
   * @brief Apply a block preconditioner to a given vector
   *
   * @param P Pointer to used preconditioner
   * @param r Vector to apply the preconditioner on
   */
  void apply_preconditioner(pavector r);
};

/* Matrix vector multiplicatio with the (approximate) schur complement for lu accuracy */
void mvm_schurcomplement(field alpha, bool trans, Block_HLU_Prcd *P, pavector x, pavector y);

void mvm_schurcomplement_rbffd(field alpha, bool trans, Block_HLU_Prcd *P, pavector x, pavector y);

real norm2diff_lr_schurcomplement_hmatrix(Block_HLU_Prcd *P, pchmatrix LR);

real norm2diff_lr_schurcomplement_hmatrix_rbffd(Block_HLU_Prcd *P, pchmatrix LR);

#endif