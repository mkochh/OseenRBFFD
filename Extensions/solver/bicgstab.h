#include "solver.h"

#ifndef BICGSTAB_HEADER
#define BICGSTAB_HEADER

class BiCGStab: public PrcdSolver
{
public: 
  uint maxiter;
  real tol;

  BiCGStab(real tol, uint maxiter);
  ~BiCGStab();

  uint solve_avector(void *A, addeval_t addeval_A, pavector b, pavector x);
  uint psolve_avector(void *A, addeval_t addeval_A, Preconditioner *P, pavector b, pavector x);
  uint psolve_avector_eigen(void *A, addeval_t addeval_A, Preconditioner *P, pavector rhs, pavector x);
};

#endif