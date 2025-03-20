#include "gmres.h"

GMRes::GMRes(real tol, uint maxiter, uint kmax)
{
  this->tol = tol;
  this->maxiter = maxiter;
  this->kmax = kmax;
}

GMRes::~GMRes()
{
  this->tol = 0.0;
  this->maxiter = 0;
  this->kmax = 0;
}

uint GMRes::solve_avector(void *A, addeval_t addeval_A, pavector b, pavector x)
{
  return solve_gmres_avector(A, addeval_A, b, x, this->tol, this->maxiter, this->kmax);
}

uint GMRes::psolve_avector(void *A, addeval_t addeval_A, Preconditioner *P, pavector b, pavector x)
{
  return solve_pgmres_avector(A, addeval_A, (prcd_t) P->prcd_function, P, b, x, this->tol, this->maxiter, this->kmax);
}