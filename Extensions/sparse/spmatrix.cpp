/**
 * \file spmatrix.cpp
 * \author Jonas Grams (jonas.grams@tuhh.de)
 * \brief
 * \version 0.1
 * \date 12-05-2022
 *
 */

#include <stdio.h>

#include "spmatrix.h"

pspmatrix new_spmatrix(psparsematrix A, psparsematrix *B, psparsematrix C)
{

  pspmatrix sp;

  sp = (pspmatrix)allocmem(sizeof(spmatrix));

  sp->A = A;
  sp->B = B;
  sp->C = C;

  return sp;
}

phspmatrix new_hspmatrix(phmatrix A, phmatrix *B, phmatrix C)
{

  phspmatrix hsp;

  hsp = (phspmatrix)allocmem(sizeof(hspmatrix));

  hsp->A = A;
  hsp->B = B;
  hsp->C = C;

  return hsp;
}

void del_spmatrix(pspmatrix sp)
{

  uint i;

  assert(sp != NULL);
  if (sp->A)
    del_sparsematrix(sp->A);
  if (sp->B)
  {
    for (i = 0; i < 3; i++)
      del_sparsematrix(sp->B[i]);
    freemem(sp->B);
  }
  if (sp->C)
    del_sparsematrix(sp->C);

  freemem(sp);
}

void del_hspmatrix(phspmatrix sp)
{

  uint i;

  assert(sp != NULL);
  if (sp->A)
    del_hmatrix(sp->A);
  if (sp->B)
  {
    for (i = 0; i < 3; i++)
      del_hmatrix(sp->B[i]);
    freemem(sp->B);
  }
  if (sp->C)
    del_hmatrix(sp->C);

  freemem(sp);
}

void del_spmatrix2(pspmatrix sp)
{

  uint i;

  assert(sp != NULL);
  if (sp->A)
    del_sparsematrix(sp->A);
  if (sp->B)
  {
    for (i = 0; i < 6; i++)
      del_sparsematrix(sp->B[i]);
    freemem(sp->B);
  }
  if (sp->C)
    del_sparsematrix(sp->C);

  freemem(sp);
}

void del_hspmatrix2(phspmatrix sp)
{

  uint i;

  assert(sp != NULL);
  if (sp->A)
    del_hmatrix(sp->A);
  if (sp->B)
  {
    for (i = 0; i < 6; i++)
      del_hmatrix(sp->B[i]);
    freemem(sp->B);
  }
  if (sp->C)
    del_hmatrix(sp->C);

  freemem(sp);
}

void addeval_spmatrix_avector(field alpha, pspmatrix sp, pcavector x, pavector y)
{

  pavector y1[3], y2;
  pavector x1[3], x2;
  uint rows_A, rows_B, cols_A;
  uint d;

  rows_A = sp->A->rows;
  rows_B = sp->B[0]->rows;
  cols_A = sp->A->cols;

  assert(x->dim == 3 * cols_A + rows_B);
  assert(y->dim == 3 * rows_A + rows_B);

  for (d = 0; d < 3; d++)
    x1[d] = new_sub_avector((pavector)x, cols_A, d * cols_A);
  x2 = new_sub_avector((pavector)x, rows_B, 3 * cols_A);
  for (d = 0; d < 3; d++)
    y1[d] = new_sub_avector(y, rows_A, d * rows_A);
  y2 = new_sub_avector(y, rows_B, 3 * rows_A);

  if (sp->A)
    for (d = 0; d < 3; d++)
      mvm_sparsematrix_avector(alpha, false, sp->A, x1[d], y1[d]);
  if (sp->B)
  {
    for (d = 0; d < 3; d++)
    {
      mvm_sparsematrix_avector(alpha, true, sp->B[d], x2, y1[d]);
      mvm_sparsematrix_avector(alpha, false, sp->B[d], x1[d], y2);
    }
  }
  if (sp->C)
    mvm_sparsematrix_avector(-1.0 * alpha, false, sp->C, x2, y2);

  for (d = 0; d < 3; d++)
    del_avector((pavector)x1[d]);
  del_avector((pavector)x2);
  for (d = 0; d < 3; d++)
    del_avector(y1[d]);
  del_avector(y2);
}

void addeval_spmatrix_avector_rbffd(field alpha, pspmatrix sp, pcavector x, pavector y)
{

  pavector y1[3], y2;
  pavector x1[3], x2;
  uint rows_A, rows_B, cols_A;
  uint d;

  rows_A = sp->A->rows;
  rows_B = sp->B[0]->rows;
  cols_A = sp->A->cols;

  assert(x->dim == 3 * cols_A + rows_B);
  assert(y->dim == 3 * rows_A + rows_B);

  for (d = 0; d < 3; d++)
    x1[d] = new_sub_avector((pavector)x, cols_A, d * cols_A);
  x2 = new_sub_avector((pavector)x, rows_B, 3 * cols_A);
  for (d = 0; d < 3; d++)
    y1[d] = new_sub_avector(y, rows_A, d * rows_A);
  y2 = new_sub_avector(y, rows_B, 3 * rows_A);

  if (sp->A)
    for (d = 0; d < 3; d++)
      mvm_sparsematrix_avector(alpha, false, sp->A, x1[d], y1[d]);
  if (sp->B)
  {
    for (d = 0; d < 3; d++)
    {
      mvm_sparsematrix_avector(alpha, true, sp->B[d+3], x2, y1[d]);
      mvm_sparsematrix_avector(alpha, false, sp->B[d], x1[d], y2);
    }
  }
  if (sp->C)
    mvm_sparsematrix_avector(alpha, false, sp->C, x2, y2);
    // mvm_sparsematrix_avector(-1.0 * alpha, false, sp->C, x2, y2);

  for (d = 0; d < 3; d++)
    del_avector((pavector)x1[d]);
  del_avector((pavector)x2);
  for (d = 0; d < 3; d++)
    del_avector(y1[d]);
  del_avector(y2);
}

void addeval_hspmatrix_avector(field alpha, phspmatrix sp, pcavector x, pavector y)
{

  pavector y1[3], y2;
  pavector x1[3], x2;
  uint rows_A, rows_B, cols_A;
  uint d;

  rows_A = sp->A->rc->size;
  rows_B = sp->B[0]->rc->size;
  cols_A = sp->A->cc->size;

  assert(x->dim == 3 * cols_A + rows_B);
  assert(y->dim == 3 * rows_A + rows_B);

  for (d = 0; d < 3; d++)
    x1[d] = new_sub_avector((pavector)x, cols_A, d * cols_A);
  x2 = new_sub_avector((pavector)x, rows_B, 3 * cols_A);
  for (d = 0; d < 3; d++)
    y1[d] = new_sub_avector(y, rows_A, d * rows_A);
  y2 = new_sub_avector(y, rows_B, 3 * rows_A);

  if (sp->A)
    for (d = 0; d < 3; d++)
      addeval_hmatrix_avector(alpha, sp->A, x1[d], y1[d]);
  if (sp->B)
  {

    for (d = 0; d < 3; d++)
    {

      addeval_hmatrix_avector(alpha, sp->B[d], x2, y1[d]);
      addevaltrans_hmatrix_avector(alpha, sp->B[d], x1[d], y2);
    }
  }
  if (sp->C)
    addeval_hmatrix_avector(-1.0 * alpha, sp->C, x2, y2);

  for (d = 0; d < 3; d++)
    del_avector(x1[d]);
  del_avector(x2);
  for (d = 0; d < 3; d++)
    del_avector(y1[d]);
  del_avector(y2);
}

void addevaltrans_spmatrix_avector(field alpha, pspmatrix sp, pcavector x, pavector y)
{

  pavector y1[3], y2;
  pavector x1[3], x2;
  uint rows_A, rows_B, cols_A;
  uint d;

  rows_A = sp->A->rows;
  rows_B = sp->B[0]->rows;
  cols_A = sp->A->cols;

  for (d = 0; d < 3; d++)
    x1[d] = new_sub_avector((pavector)x, cols_A, d * cols_A);
  x2 = new_sub_avector((pavector)x, rows_B, 3 * cols_A);
  for (d = 0; d < 3; d++)
    y1[d] = new_sub_avector(y, rows_A, d * rows_A);
  y2 = new_sub_avector(y, rows_B, 3 * rows_A);

  if (sp->A)
    for (d = 0; d < 3; d++)
      mvm_sparsematrix_avector(alpha, true, sp->A, x1[d], y1[d]);
  if (sp->B)
  {
    for (d = 0; d < 3; d++)
    {

      mvm_sparsematrix_avector(alpha, true, sp->B[d], x2, y1[d]);
      mvm_sparsematrix_avector(alpha, false, sp->B[d], x1[d], y2);
    }
  }
  if (sp->C)
    mvm_sparsematrix_avector(-1.0 * alpha, true, sp->C, x2, y2);

  for (d = 0; d < 3; d++)
    del_avector(x1[d]);
  del_avector(x2);
  for (d = 0; d < 3; d++)
    del_avector(y1[d]);
  del_avector(y2);
}

void addevaltrans_hspmatrix_avector(field alpha, phspmatrix sp, pcavector x, pavector y)
{

  pavector y1[3], y2;
  pavector x1[3], x2;
  uint rows_A, rows_B, cols_A;
  uint d;

  rows_A = sp->A->rc->size;
  rows_B = sp->B[0]->rc->size;
  cols_A = sp->A->cc->size;

  for (d = 0; d < 3; d++)
    x1[d] = new_sub_avector((pavector)x, cols_A, d * cols_A);
  x2 = new_sub_avector((pavector)x, rows_B, 3 * cols_A);
  for (d = 0; d < 3; d++)
    y1[d] = new_sub_avector(y, rows_A, d * rows_A);
  y2 = new_sub_avector(y, rows_B, 3 * rows_A);

  if (sp->A)
    for (d = 0; d < 3; d++)
      addevaltrans_hmatrix_avector(alpha, sp->A, x1[d], y1[d]);
  if (sp->B)
  {
    for (d = 0; d < 3; d++)
    {

      addevaltrans_hmatrix_avector(alpha, sp->B[d], x2, y1[d]);
      addeval_hmatrix_avector(alpha, sp->B[d], x1[d], y2);
    }
  }
  if (sp->C)
    addevaltrans_hmatrix_avector(-1.0 * alpha, sp->C, x2, y2);

  for (d = 0; d < 3; d++)
    del_avector(x1[d]);
  del_avector(x2);
  for (d = 0; d < 3; d++)
    del_avector(y1[d]);
  del_avector(y2);
}

void mvm_spmatrix_avector(field alpha, bool trans, pspmatrix sp, pcavector x, pavector y)
{

  if (trans)
    addevaltrans_spmatrix_avector(alpha, sp, x, y);
  else
    addeval_spmatrix_avector(alpha, sp, x, y);
}

void mvm_hspmatrix_avector(field alpha, bool trans, phspmatrix sp, pcavector x, pavector y)
{

  if (trans)
    addevaltrans_hspmatrix_avector(alpha, sp, x, y);
  else
    addeval_hspmatrix_avector(alpha, sp, x, y);
}
