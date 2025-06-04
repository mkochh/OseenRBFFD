#include "hlu.h"
#include "../cluster/oseencluster.h"
#include "../auxiliaries/aux_medusa.hpp"
#include "../auxiliaries/aux_h2lib.h"

HLU_Options::HLU_Options(pspmatrix K, uint clfu, uint clfp,
                         clustermode modeu, clustermode modep, 
                         void *eta_vel, void* eta_grad, void*eta_schur,
                         admissible admu, admissible admpu, admissible admp, ptruncmode tm,
                         real eps_vel, real eps_mul, real eps_schur)
{

  this->K = K;

  this->clfu = clfu;
  this->clfp = clfp;

  this->modeu = modeu;
  this->modep = modep;

  this->eta_vel = eta_vel;
  this->eta_grad = eta_grad;
  this->eta_schur = eta_schur;

  this->admu = admu;
  this->admpu = admpu;
  this->admp = admp;

  this->tm = tm;

  this->eps_vel = eps_vel;
  this->eps_mul = eps_mul;
  this->eps_schur = eps_schur;
}

Block_HLU_Prcd::Block_HLU_Prcd(mm::DomainDiscretization<mm::Vec3d> d_u,
                                mm::DomainDiscretization<mm::Vec3d> d_div,
                                mm::Range<int> p_idxs_in_d_div,
                                HLU_Options *opt,
                                Block_HLU_Times *times,
                                psparsematrix symA,
                                BlockPrcdType prcd_type)
{
  phmatrix Ahat, Bhat1, Bhat2, Shat;
  pclustergeometry cgv, cgp;
  pcluster rootv, rootp, rootpconstraint, rootconstraint, rootrootv;
  pblock blocka, blockb, blocks;
  uint *idxv, *idxp, *flag, *idxpconstraint;
  pstopwatch sw = new_stopwatch();

  flag = new uint[d_u.all().size()];

  for (uint i = 0; i < static_cast<uint>(d_u.all().size()); i++)
    flag[i] = 0;

  mm::Range<int> p_idxs_in_d_p = p_idxs_in_d_div;
  int Nu = p_idxs_in_d_div[0];
  for (int i = 0; i < p_idxs_in_d_div.size(); i++)
    p_idxs_in_d_p[i] -= Nu;
  
  idxv = idxMM2H(d_u.all());
  idxp = idxMM2H(p_idxs_in_d_p);

  assert(p_idxs_in_d_p.size() == p_idxs_in_d_div.size());
  
  cgv = build_clustergeometry_medusa(d_u, d_u.all());
  cgp = build_clustergeometry_medusa(d_div, p_idxs_in_d_div);

  start_stopwatch(sw);
  rootp = build_cluster(cgp, p_idxs_in_d_p.size(), idxp, opt->clfp, opt->modep);
  times->pressure_cluster = stop_stopwatch(sw);

  start_stopwatch(sw);
  rootv = build_adaptive_dd_cluster(cgv, d_u.all().size(), idxv, opt->clfu, symA, 3, flag);
  times->velocity_cluster = stop_stopwatch(sw);

  delete[] flag;
  
  // construct blocks such that pressure constraint is included
  rootrootv = new_cluster(rootv->size, idxv, 1, 3); // auxiliary cluster to make construction of blockb easier
  rootrootv->son[0] = rootv;
  blocka = new_block(rootrootv, rootrootv, false, 1, 1);
  blocka->son[0] = build_nonstrict_block(rootv, rootv, opt->eta_vel, opt->admu);
  update_block(blocka);
  
  idxpconstraint = (uint*)malloc((p_idxs_in_d_p.size()+1)*sizeof(uint)); // p_idxs.size() = Anzahl der pressure nodes
  for (int i = 0; i < p_idxs_in_d_p.size(); i++) { idxpconstraint[i] = idxp[i]; }
  idxpconstraint[p_idxs_in_d_p.size()] = p_idxs_in_d_p.size(); // last index for pressure constraint
  
  rootconstraint = new_cluster(1, idxpconstraint+(p_idxs_in_d_p.size()), 0, 3); // cluster containing only the pressure constraint
  rootpconstraint = new_cluster(p_idxs_in_d_p.size()+1, idxpconstraint, 2, 3); // cluster containing pressure nodes and constraint
  rootpconstraint->son[0] = rootp;
  rootpconstraint->son[1] = rootconstraint;
  blocks = new_block(rootpconstraint, rootpconstraint, false, 2, 2);
  blocks->son[0] = build_nonstrict_block(rootp, rootp, opt->eta_schur, opt->admp);
  blocks->son[1] = new_block(rootconstraint, rootp, false, 0, 0);
  blocks->son[2] = new_block(rootp, rootconstraint, false, 0, 0);
  blocks->son[3] = new_block(rootconstraint, rootconstraint, false, 0, 0);
  update_block(blocks);

  blockb = new_block(rootpconstraint, rootrootv, false, 2, 1);
  blockb->son[0] = build_nonstrict_block(rootp, rootv, opt->eta_grad, opt->admpu);
  blockb->son[1] = new_block(rootconstraint, rootv, false, 0, 0);
  update_block(blockb);

  std::cout << "block structure for H-matrices build" << std::endl;

  times->velocity_compression = 0.0;
  Ahat = build_from_block_hmatrix(blocka, 0);
  copy_sparsematrix_to_hmatrix(opt->K->A, Ahat);

  Shat = build_from_block_hmatrix(blocks, 0);
  copy_sparsematrix_to_hmatrix(opt->K->C, Shat);

  start_stopwatch(sw);
  lrdecomp_hmatrix(Ahat, opt->tm, opt->eps_vel);
  times->velocity_lu = stop_stopwatch(sw);

  std::cout << "sparse matrices copied to H-matrices" << std::endl;
    
  times->grad_lower_solve = 0.0;
  times->grad_upper_solve = 0.0;
  times->grad_schur_multiplication = 0.0;
  times->gradient_compression = 0.0;
  #ifdef USE_OPENMP
  #pragma omp parallel for shared(Shat) private(Bhat1, Bhat2)
  #endif
  for (uint d = 0; d < 3; d++)
  {
    std::cout << "computing schur complement summands" << std::endl;
    #ifdef USE_OPENMP
    #pragma omp parallel for
    #endif
    for (int i = 0; i < 2; i++) {
      if (i == 0) {
        Bhat1 = build_from_block_hmatrix(blockb, 0);
        copy_sparsematrix_to_hmatrix(opt->K->B[d], Bhat1);

        start_stopwatch(sw);
        triangularinvmul_hmatrix(false, false, true, Ahat, opt->tm, opt->eps_mul, true, Bhat1);
        times->grad_upper_solve += stop_stopwatch(sw);
      } else {
        Bhat2 = build_from_block_hmatrix(blockb, 0);
        copy_sparsematrix_to_hmatrix(opt->K->B[d+3], Bhat2);

        start_stopwatch(sw);
        triangularinvmul_hmatrix(true, true, false, Ahat, opt->tm, opt->eps_mul, true, Bhat2);
        times->grad_lower_solve += stop_stopwatch(sw);
      }
    }

    start_stopwatch(sw);
    phmatrix Temp;
    Temp = build_from_block_hmatrix(blocks, 0);
    clear_hmatrix(Temp);
    addmul_hmatrix(-1.0, false, Bhat1, true, Bhat2, opt->tm, opt->eps_mul, Temp);

    std::cout << "adding schur complement summands" << std::endl;
    #ifdef USE_OPENMP
    #pragma omp critical
    #endif
    {
      add_hmatrix(1.0, Temp, opt->tm, opt->eps_mul, Shat);
    }
    times->grad_schur_multiplication += stop_stopwatch(sw);

    del_hmatrix(Bhat1);
    del_hmatrix(Bhat2);
    del_hmatrix(Temp);
  }

  // Set total times for the schur complement computation
  times->schur_computation = times->grad_schur_multiplication + times->grad_lower_solve + times->grad_upper_solve;

  start_stopwatch(sw);
  lrdecomp_hmatrix(Shat, opt->tm, opt->eps_schur);
  times->schur_lu = stop_stopwatch(sw);

  this->A = Ahat;
  this->S = Shat;
  this->B = opt->K->B;
  this->C = opt->K->C;
  this->block_velocity = blocka;
  this->block_grad = blockb;
  this->block_pressure = blocks;
  this->idxu = idxv;
  this->idxp = idxpconstraint;
  this->root_velocity = rootrootv;
  this->root_pressure = rootpconstraint;
  this->prcd_type = prcd_type;

  del_stopwatch(sw);
  del_clustergeometry(cgv);
  del_clustergeometry(cgp);
}

Block_HLU_Prcd::~Block_HLU_Prcd()
{

  del_block(this->block_velocity);
  del_block(this->block_pressure);
  del_block(this->block_grad);

  del_cluster(this->root_velocity);
  delete[] this->root_pressure->son[0]->idx;
  del_cluster(this->root_pressure);

  delete[] this->idxu;
  delete[] this->idxp;

  del_hmatrix(this->A);
  del_hmatrix(this->S);
}

void Block_HLU_Prcd::apply_preconditioner(pavector r)
{

  pavector *r1, r2;
  uint n, m, d;

  r1 = new pavector[3];

  n = this->A->rc->size;
  m = this->S->rc->size;

  assert(3 * n + m == r->dim);

  for (d = 0; d < 3; d++)
    r1[d] = new_sub_avector(r, n, d * n);
  r2 = new_sub_avector(r, m, 3 * n);

  /* r1 <- A^{-1}*r1 */
  for (d = 0; d < 3; d++)
    lrsolve_hmatrix_avector(false, this->A, r1[d]);

  // Apply lower triangular part
  if (this->prcd_type == BLOCK_TRIANGULAR)
  {
    for(d = 0; d < 3; d++)
      addeval_sparsematrix_avector(-1.0, this->B[d], r1[d], r2);
  }

  /* r2 <- -S^{-1}*r2 */
  lrsolve_hmatrix_avector(false, this->S, r2);
  // scale_avector(-1.0, r2); // why?

  for (d = 0; d < 3; d++)
    del_avector(r1[d]);

  delete[] r1;

  del_avector(r2);
}

void mvm_schurcomplement(field alpha, bool trans, Block_HLU_Prcd *P, pavector x, pavector y)
{
  pavector tmp;
  tmp = new_zero_avector(P->A->rc->size);

  for(uint d = 0; d < 3; d++)
  {
    clear_avector(tmp);
    mvm_sparsematrix_avector(1.0, true, P->B[d], x, tmp);
    lrsolve_hmatrix_avector(trans, P->A, tmp);
    mvm_sparsematrix_avector(-alpha, false, P->B[d], tmp, y);
  }

  del_avector(tmp);
}

void mvm_schurcomplement_rbffd(field alpha, bool trans, Block_HLU_Prcd *P, pavector x, pavector y)
{
  pavector tmp;
  tmp = new_zero_avector(P->A->rc->size);

  if (trans)
  {
    for(uint d = 0; d < 3; d++)
    {
      clear_avector(tmp);
      mvm_sparsematrix_avector(1.0, true, P->B[d], x, tmp);
      lrsolve_hmatrix_avector(trans, P->A, tmp);
      mvm_sparsematrix_avector(-alpha, false, P->B[d+3], tmp, y);
    }
    mvm_sparsematrix_avector(1.0, true, P->C, x, y);
  } else {
    for(uint d = 0; d < 3; d++)
    {
      clear_avector(tmp);
      mvm_sparsematrix_avector(1.0, true, P->B[d+3], x, tmp);
      lrsolve_hmatrix_avector(trans, P->A, tmp);
      mvm_sparsematrix_avector(-alpha, false, P->B[d], tmp, y);
    }
    mvm_sparsematrix_avector(1.0, false, P->C, x, y);
  }
  del_avector(tmp);
}

real
norm2diff_lr_schurcomplement_hmatrix(Block_HLU_Prcd *P, pchmatrix LR)
{
  return norm2diff_pre_matrix((mvm_t) mvm_schurcomplement, (void *) P,
			      (prcd_t) lreval_n_hmatrix_avector,
			      (prcd_t) lreval_t_hmatrix_avector, (void *) LR,
			      P->S->rc->size, P->S->cc->size);
}

real
norm2diff_lr_schurcomplement_hmatrix_rbffd(Block_HLU_Prcd *P, pchmatrix LR)
{
  return norm2diff_pre_matrix((mvm_t) mvm_schurcomplement_rbffd, (void *) P,
			      (prcd_t) lreval_n_hmatrix_avector,
			      (prcd_t) lreval_t_hmatrix_avector, (void *) LR,
			      P->S->rc->size, P->S->cc->size);
}