#include "../../H2Lib/h2lib.h"

#ifndef ADMISSIBLE_HEADER
#define ADMISSIBLE_HEADER

typedef struct
{
  uint nmin;
  psparsematrix sp;
  
} adm_sparse_data;

/* ************************************
 * Coupled admissibility conditions   *
 **************************************/

bool
admissible_coupled_cluster(pcluster s, pcluster t, void *data);

bool
admissible_coupled_cluster_v2(pcluster rc, pcluster cc, void *data);

/****************************************
 * Weaker admissibility conditions      *
 ****************************************/

bool
admissible_hodlr(pcluster s, pcluster t, void *data);

bool 
admissible_sparse(pcluster s, pcluster t, void *data);

bool
admissible_weak(pcluster s, pcluster t, void *data);

bool
admissible_weak_rbffd(pcluster s, pcluster t, void *data);

/* *****************************************
 * DD with weaker admissibility conditions *
 *******************************************/

bool
admissible_dd_sparse(pcluster s, pcluster t, void *data);

bool
admissible_dd_weak(pcluster s, pcluster t, void *data);

bool
admissible_dd_weak_rbffd(pcluster s, pcluster t, void *data);

/* **********************************************
 * Coupled with weaker admissibility conditions *
 ************************************************/

bool
admissible_coupled_sparse(pcluster s, pcluster t, void *data);

bool
admissible_coupled_weak(pcluster s, pcluster t, void *data);

bool
admissible_coupled_hodlr(pcluster s, pcluster t, void *data);

/* *****************************************
 * IA admissibility conditions *
 *******************************************/

bool
admissible_ia_cluster(pcluster s, pcluster t, void *data);

bool
admissible_ia_sparse(pcluster s, pcluster t, void *data);

bool
admissible_ia_weak(pcluster s, pcluster t, void *data);

/* *****************************************
 * IA coupled admissibility conditions *
 *******************************************/

bool
admissible_coupled_ia_cluster(pcluster s, pcluster t, void *data);

bool
admissible_coupled_ia_sparse(pcluster s, pcluster t, void *data);

bool
admissible_coupled_ia_weak(pcluster s, pcluster t, void *data);

#endif