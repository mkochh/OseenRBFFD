/**
 * \file spmatrix.h
 * \author Jonas Grams (jonas.grams@tuhh.de)
 * \brief Representation of a saddle point matrix.
 * \version 1.0
 * \date 12-05-2022
 *
 */

#include "../../H2Lib/h2lib.h"

#ifndef SPMATRIX_HEADER
#define SPMATRIX_HEADER

/* Representation of a saddle point matrix, i.e. of a 2x2 block matrix
   with (1,1) block A, (2,1) block B, (2,1) block B^T, (2,2) block 0 */
struct _spmatrix
{

  psparsematrix A;  // (1,1) block of the matrix
  psparsematrix *B; // (2,1) block of the matrix, and transposed of the (1,2) block
  psparsematrix C;  // Stabilization part
};

/* Representation of a saddle point matrix with H-Matrix blocks */
struct _hspmatrix
{

  phmatrix A;
  phmatrix *B;
  phmatrix C;
};

typedef struct _spmatrix spmatrix;
typedef spmatrix *pspmatrix;
typedef const pspmatrix pcspmatrix;

typedef struct _hspmatrix hspmatrix;
typedef hspmatrix *phspmatrix;
typedef const phspmatrix pchspmatrix;

/* Constructor and destructor for a saddale point matrix */
pspmatrix new_spmatrix(psparsematrix A, psparsematrix *B, psparsematrix C);
phspmatrix new_hspmatrix(phmatrix A, phmatrix *B, phmatrix C);

void del_spmatrix(pspmatrix);
void del_hspmatrix(phspmatrix);

void del_spmatrix2(pspmatrix);
void del_hspmatrix2(phspmatrix);

/* Matrix vector multiplication for saddle point matrices */
void addeval_spmatrix_avector(field alpha, pspmatrix sp, pcavector x, pavector y);
void addeval_spmatrix_avector_rbffd(field alpha, pspmatrix sp, pcavector x, pavector y);
void addeval_hspmatrix_avector(field alpha, phspmatrix sp, pcavector x, pavector y);

void addevaltrans_spmatrix_avector(field alpha, pspmatrix sp, pcavector x, pavector y);
void addevaltrans_hspmatrix_avector(field alpha, phspmatrix sp, pcavector x, pavector y);

void mvm_spmatrix_avector(field alpha, bool trans, pspmatrix sp, pcavector x, pavector y);
void mvm_hspmatrix_avector(field alpha, bool trans, phspmatrix sp, pcavector x, pavector y);
#endif