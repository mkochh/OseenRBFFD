
#include "H2Lib/h2lib.h"

/* Include extended modules */
#include "Extensions/preconditioning/prec.h"
#include "Extensions/preconditioning/block_diagonal.h"
#include "Extensions/preconditioning/hlu.h"

#include "Extensions/sparse/spmatrix.h"
#include "Extensions/sparse/sparse_compression.h"

#include "Extensions/cluster/oseencluster.h"
#include "Extensions/cluster/admissible.h"

#include "Extensions/solver/solver.h"
#include "Extensions/solver/bicgstab.h"
#include "Extensions/solver/gmres.h"

#include "Extensions/auxiliaries/aux_h2lib.h"
#include "Extensions/auxiliaries/aux_medusa.hpp"
#include "Extensions/auxiliaries/aux_eigen.hpp"
#include "Extensions/auxiliaries/polyhedron_integration.h"
#include "Extensions/auxiliaries/spmatrix_medusa.h"

#include "Extensions/oseen_weights_rbffd/set_weights.hpp"
#include "Extensions/oseen_weights_rbffd/rbffd_discretization.hpp"
