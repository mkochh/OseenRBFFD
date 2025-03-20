#include "../../H2Lib/h2lib.h"

#ifndef SPARSECMPRESSION_HEADER
#define SPARSECMPRESSION_HEADER

psparsematrix
reorder_sparsematrix(psparsematrix sp, uint *rowinv, uint *col);

void 
copy_sparsematrix_to_hmatrix(psparsematrix src, phmatrix trg);
#endif