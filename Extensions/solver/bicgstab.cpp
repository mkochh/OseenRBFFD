#include "bicgstab.h"

BiCGStab::BiCGStab(real tol, uint maxiter)
{
  this->tol = tol;
  this->maxiter = maxiter;
}

BiCGStab::~BiCGStab()
{
  this->tol = 0.0;
  this->maxiter = 0;
}

uint BiCGStab::solve_avector(void *A, addeval_t addeval_A, pavector b, pavector x)
{
  uint n, iter;
  pavector r, rt, p, a, as;
  real tol, error;

  n = x->dim;

  assert(n == b->dim);

  r = new_zero_avector(n);
  rt = new_avector(n);
  p = new_avector(n);
  a = new_avector(n);
  as = new_avector(n);

  tol = this->tol * norm2_avector(b);

  init_bicgstab(addeval_A, A, b, x, r, rt, p, a, as);
  error = norm2_avector(r);

  iter = 0;
  while (error > tol && iter + 1 != maxiter)
  {

    step_bicgstab(addeval_A, A, b, x, r, rt, p, a, as);
    error = norm2_avector(r);

    iter++;
  }

  del_avector(r);
  del_avector(rt);
  del_avector(p);
  del_avector(a);
  del_avector(as);

  return iter;
}

/* ------------------------------------------------------------
 * Preconditioned stabilized biconjugate gradient method
 * ------------------------------------------------------------ */

/* cf. Yousef Saad, Iterative Methods for Sparse Linear Systems,
 Section 7.4.2, where every multiplication with A comes along with
 a call of some preconditioning routine.
 Slight modification: the intermediate vector s is stored in r. */

static void init_pbicgstab(addeval_t addeval, void *matrix, prcd_t prcd, void *pdata, pcavector b, /* Right-hand side */
                           pavector x,                                                             /* Approximate solution */
                           pavector r,                                                             /* Residual b-Ax */
                           pavector rt,                                                            /* Adjoint residual */
                           pavector p,                                                             /* Search direction */
                           pavector a, pavector as)
{
  (void)a;
  (void)as;

  copy_avector(b, r); /* r = b - A x */
  addeval(-1.0, matrix, x, r);
  if (prcd)
    prcd(pdata, r);

  copy_avector(r, rt); /* r^* = r */

  copy_avector(r, p); /* p = r */
}

static void step_pbicgstab(addeval_t addeval, void *matrix, prcd_t prcd, void *pdata, pcavector b, /* Right-hand side */
                           pavector x,                                                             /* Approximate solution */
                           pavector r,                                                             /* Residual b-Ax */
                           pavector rt,                                                            /* Adjoint residual */
                           pavector p,                                                             /* Search direction */
                           pavector a, pavector as)
{
  field alpha, beta, omega, mu;

  (void)b;

  clear_avector(a); /* a = A p */
  addeval(1.0, matrix, p, a);
  if (prcd)
    prcd(pdata, a);

  mu = dotprod_avector(r, rt);
  alpha = mu / dotprod_avector(a, rt);

  //printf("mu = %8.7e,\talpha = %8.7e\n", mu, alpha);

  add_avector(-alpha, a, r); /* r = r - alpha a */

  clear_avector(as); /* as = A r */
  addeval(1.0, matrix, r, as);
  if (prcd)
    prcd(pdata, as);

  omega = dotprod_avector(as, r) / dotprod_avector(as, as);

  //printf("omega = %8.7e", omega);

  add_avector(alpha, p, x); /* x = x + alpha p + omega s */
  add_avector(omega, r, x);

  add_avector(-omega, as, r); /* r = r - omega as */

  beta = dotprod_avector(r, rt) / mu * alpha / omega;

  scale_avector(beta, p); /* p = r + beta (p - omega a) */
  add_avector(1.0, r, p);
  add_avector(-beta * omega, a, p);
}

uint BiCGStab::psolve_avector(void *A, addeval_t addeval_A, Preconditioner *P, pavector b, pavector x)
{
  uint n, iter;
  pavector r, rt, p, a, as;
  real tol, error;

  n = x->dim;

  assert(n == b->dim);

  r = new_zero_avector(n);
  rt = new_avector(n);
  p = new_avector(n);
  a = new_avector(n);
  as = new_avector(n);

  tol = this->tol * norm2_avector(b);

  init_pbicgstab(addeval_A, A, (prcd_t) P->prcd_function, P, b, x, r, rt, p, a, as);
  error = norm2_avector(r);

  iter = 0;
  while (error > tol && iter + 1 != maxiter)
  {
    //printf("\niter = %5d and 2-norm error of residual = %8.7e\n", iter, error);

    step_pbicgstab(addeval_A, A, (prcd_t) P->prcd_function, P, b, x, r, rt, p, a, as);
    error = norm2_avector(r);

    iter++;
  }

  del_avector(r);
  del_avector(rt);
  del_avector(p);
  del_avector(a);
  del_avector(as);

  return iter;
}

static void apply_prcd(prcd_t prcd, void* pdata, pavector a) {
  if (prcd)
    prcd(pdata, a);
}

uint BiCGStab::psolve_avector_eigen(void *A, addeval_t addeval_A, Preconditioner *P, pavector rhs, pavector x)
{
  uint n, iter;
  pavector r, r0, p, v, y, s, z, t;
  real tol2, eps2, error;

  n = x->dim;

  assert(n == rhs->dim);

  r = new_zero_avector(n);
  r0 = new_avector(n);
  p = new_avector(n);
  v = new_zero_avector(n);
  s = new_avector(n);
  y = new_avector(n);
  z = new_avector(n);
  t = new_avector(n);

  copy_avector(rhs, r); 
  addeval_A(-1.0, A, x, r); // r = rhs - A x

  copy_avector(r, p); // p = r 
  copy_avector(r, r0); // r0 = r

  real r0_sqnorm = dotprod_avector(r0,r0);
  real rhs_sqnorm = dotprod_avector(rhs,rhs);
  if(rhs_sqnorm == 0) {
    clear_avector(x);
    return 0;
  }

  real rho = 1;
  real alpha = 1;
  real w = 1;

  tol2 = this->tol * norm2_avector(rhs);
  eps2 = 1e-32;

  error = norm2_avector(r);

  iter = 0;
  while (error > tol2 && iter < maxiter)
  {
    // printf("\niter = %5d and 2-norm error of residual = %8.7e\n", iter, error);

    real rho_old = rho;

    rho = dotprod_avector(r0, r);
    if(abs(rho) < eps2*r0_sqnorm) {
      // The new residual vector became too orthogonal to the arbitrarily chosen direction r0
      // Let's restart with a new r0:
      copy_avector(rhs, r); 
      addeval_A(-1.0, A, x, r); // r = rhs - A x
      rho = r0_sqnorm = dotprod_avector(r,r);
    }
    real beta = (rho/rho_old) * (alpha/w);
    scale_avector(beta, p); 
    add_avector(1.0, r, p);
    add_avector(-beta*w, v, p); /* p = r + beta (p - omega v) */

    copy_avector(p, y);
    apply_prcd((prcd_t) P->prcd_function, P, y); // y = precond.solve(p)

    clear_avector(v);
    addeval_A(1.0, A, y, v); // v = A*y

    alpha = rho / dotprod_avector(r0,v);
    copy_avector(r,s);
    add_avector(-1.0*alpha, v, s); // s = r - alpha*v

    copy_avector(s, z);
    apply_prcd((prcd_t) P->prcd_function, P, z); // z = precond.solve(s)
    clear_avector(t);
    addeval_A(1.0, A, z, t); // t = A*z

    real tmp = dotprod_avector(t,t);
    if(tmp>real(0))
      w = dotprod_avector(t,s)/tmp;
    else
      w = real(0);

    add_avector(alpha, y, x);
    add_avector(w, z, x); // x = x + alpha * y + w * z

    copy_avector(s, r);
    add_avector(-1.0*w, t, r);

    error = norm2_avector(r);

    iter++;
  }

  del_avector(r);
  del_avector(r0);
  del_avector(p);
  del_avector(v);
  del_avector(y);
  del_avector(s);

  return iter;
}