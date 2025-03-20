#include "block_diagonal.h"

Block_Diagonal_Prcd::Block_Diagonal_Prcd(Preconditioner *prcd_A, uint n,
                                         Preconditioner *prcd_S, uint m)
{

  this->prcd_A = prcd_A;
  this->prcd_S = prcd_S;

  this->n = n;
  this->m = m;
}

Block_Diagonal_Prcd::~Block_Diagonal_Prcd()
{

}

void Block_Diagonal_Prcd::apply_preconditioner(pavector r)
{

  uint n = this->n;
  uint m = this->m;

  pavector *r1 = new pavector[3];
  pavector r2;

  for (int i = 0; i < 3; i++)
    r1[i] = new_sub_avector(r, n, i * n);
  r2 = new_sub_avector(r, m, 3 * n);

  for (int i = 0; i < 3; i++)
    this->prcd_A->apply_preconditioner(r1[i]);
  this->prcd_S->apply_preconditioner(r2);

  for (int i = 0; i < 3; i++)
    del_avector(r1[i]);

  delete[] r1;
  del_avector(r2);
}
