#include "../../H2Lib/h2lib.h"
#include "../sparse/spmatrix.h"
#include "../sparse/sparse_compression.h"
#include "../cluster/oseencluster.h"


#ifndef PREC_HEADER
#define PREC_HEADER

/**
 * @brief Base class for preconditioner.
 *        Each preconditioner should have an function implementing the
 *        application to a vector
 */
class Preconditioner
{
public:
  virtual ~Preconditioner(){};
  virtual void apply_preconditioner(pavector){};
  static void prcd_function(Preconditioner *P, pavector r)
  {
    P->apply_preconditioner(r);
  };
};

#endif