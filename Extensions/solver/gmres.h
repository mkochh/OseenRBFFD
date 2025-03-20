#include "solver.h"
#include "../../H2Lib/h2lib.h"

#ifndef GMRES_HEADER
#define GMRES_HEADER

class GMRes: public PrcdSolver
{
public: 
  uint maxiter;
  real tol;
  uint kmax;

  GMRes(real tol, uint maxiter, uint kmax);
  ~GMRes();

  uint solve_avector(void *A, addeval_t addeval_A, pavector b, pavector x);
  uint psolve_avector(void *A, addeval_t addeval_A, Preconditioner *P, pavector b, pavector x);
};

#endif