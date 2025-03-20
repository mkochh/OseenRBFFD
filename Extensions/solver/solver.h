#include "../../H2Lib/h2lib.h"
#include "../preconditioning/prec.h"

#ifndef SOLVER_HEADER
#define SOLVER_HEADER
class Solver
{
public:
  Solver(){};
  virtual ~Solver(){};

  virtual uint solve_avector(void *, addeval_t, pavector, pavector){return 0;};
};

class PrcdSolver: public Solver
{
public:
  PrcdSolver(){};
  virtual ~PrcdSolver(){};

  virtual uint solve_avector(void *, addeval_t, pavector, pavector){return 0;};
  virtual uint psolve_avector(void *, addeval_t, Preconditioner *, pavector, pavector){return 0;};
};

#endif