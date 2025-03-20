#include "admissible.h"
#include "../preconditioning/hlu.h"

/* ************************************
 * Coupled admissibility conditions   *
 **************************************/

bool admissible_coupled_cluster(pcluster s, pcluster t, void *data)
{
  bool b, a;

  a = admissible_2_min_cluster(s, t, data);

  if (a == true)
    b = true;
  else
  {
    a = (s->associated == t);
    if ((t->type == 1 && s->associated->type == t->type && a == false /*&& s != t */))
      b = true;
    else
      b = false;
  }
  return b;
}

bool admissible_coupled_cluster_v2(pcluster rc, pcluster cc, void *data)
{
  bool a = admissible_2_cluster(rc, cc, data);

  if (a == true)
  {
    return true;
  }
  else if(rc->type == 1 && cc->type == 1 && rc != cc)
  {
    return true;
  }
  else if(rc->type == 4 && cc->type == 1)
  {
    return true;
  }
  else if(rc->type == 1 && cc->type == 4)
  {
    return true;
  }
  else if(rc->type == 3 && cc->type == 1 && rc->associated != cc && rc->associated != NULL)
  {
    return true;
  }
  else if(rc->type == 1 && cc->type == 3 && cc->associated != rc && cc->associated != NULL)
  {
    return true;
  }
  else if(rc->type == 3 && cc->type == 3 && rc != cc)
  {
    return true;
  }
  else
  {
    return false;
  }
}

/****************************************
 * Weaker admissibility conditions      *
 ****************************************/

bool admissible_hodlr(pcluster s, pcluster t, void *data)
{
  (void) data;

  if(s == t)
  {
    return false;
  }

  return true;
}

bool admissible_sparse(pcluster s, pcluster t, void *data)
{
  adm_sparse_data * adm_data = *(adm_sparse_data **) data;
  uint nmin = adm_data->nmin;
  psparsematrix sp = adm_data->sp;

  if(s == t)
  {
    return false;
  }
  
  uint *flag = new uint[sp->cols];
  for(uint i = 0; i < sp->cols; i++)
  {
    flag[i] = 0;
  }

  for(uint i = 0; i < s->size; i++)
  {
    uint ii = s->idx[i];
    for(uint k = sp->row[ii]; k < sp->row[ii+1]; k++)
    {
      uint j = sp->col[k];
      flag[j] = 1;
    }
  }

  uint counter = 0;
  for(uint i = 0; i < t->size; i++)
  {
    uint ii = t->idx[i];
    if(flag[ii] == 1)
    {
      counter++;
    }
  }

  delete[] flag;

  if(counter < nmin)
  {
    return true;
  }

  return false;
}

bool admissible_weak(pcluster s, pcluster t, void *data)
{
  real dist;

  (void) data;

  dist = getdist_2_cluster(s, t);

  return (dist > 0);
}

bool admissible_weak_rbffd(pcluster s, pcluster t, void *data)
{
  real dist;
  real rad = *(real *) data;

  dist = getdist_2_cluster(s, t);

  return (dist > rad);
}

/* *****************************************
 * DD with weaker admissibility conditions *
 *******************************************/

bool admissible_dd_sparse(pcluster s, pcluster t, void *data)
{
  bool a;

  a = admissible_sparse(s, t, data);

  if(a == true)
  {
    return true;
  }
  else{
    a = (s == t);
    if((s->type == 1) && (s->type == t->type) && (a == false))
    {
      return true;
    }
    else
    {
      return false;
    }
  }
}

bool admissible_dd_weak(pcluster s, pcluster t, void *data)
{
  bool a;

  a = admissible_weak(s, t, data);

  if(a == true)
  {
    return true;
  }
  else{
    a = (s == t);
    if((s->type == 1) && (s->type == t->type) && (a == false))
    {
      return true;
    }
    else
    {
      return false;
    }
  }
}

bool admissible_dd_weak_rbffd(pcluster s, pcluster t, void *data)
{
  bool a;

  a = admissible_weak_rbffd(s, t, data);

  if(a == true)
  {
    return true;
  }
  else{
    a = (s == t);
    if((s->type == 1) && (s->type == t->type) && (a == false))
    {
      return true;
    }
    else
    {
      return false;
    }
  }
}

/* **********************************************
 * Coupled with weaker admissibility conditions *
 ************************************************/

bool admissible_coupled_sparse(pcluster s, pcluster t, void *data)
{
  bool b, a;

  assert(t->type > 0);

  a = admissible_sparse(s, t, data);

  if (a == true)
    b = true;
  else
  {
    a = (s->associated == t);
    b = (s == t->associated);

    assert(a == b);
    if ((t->type == 1 && s->associated->type == t->type && a == false))
      b = true;
    else
      b = false;
  }
  return b;
}

bool admissible_coupled_weak(pcluster s, pcluster t, void *data)
{
  bool b, a;

  assert(t->type > 0);

  (void) data;

  a = admissible_weak(s, t, 0);

  if (a == true)
    b = true;
  else
  {
    a = (s->associated == t);
    b = (s == t->associated);

    assert(a == b);
    if ((t->type == 1 && s->associated->type == t->type && a == false))
      b = true;
    else
      b = false;
  }
  return b;
}

bool admissible_coupled_hodlr(pcluster s, pcluster t, void *data)
{
  bool b, a;

  assert(t->type > 0);

  (void) data;
  a = admissible_hodlr(s->associated, t, 0);
  b = admissible_hodlr(s, t->associated, 0);

  if ((a == true) || (b == true))
    return true;
  else
  {
    a = (s->associated == t);
    b = (s == t->associated);

    assert(a == b);
    if ((t->type == 1 && s->associated->type == t->type && a == false))
      return true;
    else
      return false;
  }
}

/* *****************************************
 * IA admissibility conditions *
 *******************************************/

bool admissible_ia_cluster(pcluster s, pcluster t, void *data)
{
    bool a = admissible_dd_cluster(s, t, data);
    bool b;

    if (a == true)
        b = true;
    else if (s->type == 2 && t->type == 1)
    {
        a = (s->left == t) ||
            (s->right == t) ||
            (s == t);

        if (a == false)
            b = true;
        else
            b = false;
    }
    else if(t->type == 2 && s->type == 1)
    {
        a = (t->left == s) ||
            (t->right == s) ||
            (s == t);

        b = a ? false : true;
    }
    else
    {
        b = false;
    }

    return b;
}

bool admissible_ia_sparse(pcluster s, pcluster t, void *data)
{
    bool a = admissible_dd_sparse(s, t, data);
    bool b;

    if (a == true)
        b = true;
    else if (s->type == 2 && t->type == 1)
    {
        a = (s->left == t) ||
            (s->right == t) ||
            (s == t);

        if (a == false)
            b = true;
        else
            b = false;
    }
    else if(t->type == 2 && s->type == 1)
    {
        a = (t->left == s) ||
            (t->right == s) ||
            (s == t);

        b = a ? false : true;
    }
    else
    {
        b = false;
    }

    return b;
}

bool admissible_ia_weak(pcluster s, pcluster t, void *data)
{
    bool a = admissible_dd_weak(s, t, data);
    bool b;

    if (a == true)
        b = true;
    else if (s->type == 2 && t->type == 1)
    {
        a = (s->left == t) ||
            (s->right == t) ||
            (s == t);

        if (a == false)
            b = true;
        else
            b = false;
    }
    else if(t->type == 2 && s->type == 1)
    {
        a = (t->left == s) ||
            (t->right == s) ||
            (s == t);

        b = a ? false : true;
    }
    else
    {
        b = false;
    }

    return b;
}

/* *****************************************
 * IA coupled admissibility conditions *
 *******************************************/

bool admissible_coupled_ia_cluster(pcluster s, pcluster t, void *data)
{
    bool b, a;

    a = admissible_2_min_cluster(s, t, data);

    if (a == true)
        b = true;
    else
    {
        b = admissible_ia_cluster(s->associated, t, data);
    }
    return b;
}

bool admissible_coupled_ia_sparse(pcluster s, pcluster t, void *data)
{
    bool b, a;

    a = admissible_sparse(s, t, data);

    if (a == true)
        b = true;
    else
    {
        b = admissible_ia_sparse(s->associated, t, data);
    }
    return b;
}

bool admissible_coupled_ia_weak(pcluster s, pcluster t, void *data)
{
    bool b, a;

    a = admissible_weak(s, t, data);

    if (a == true)
        b = true;
    else
    {
        b = admissible_ia_weak(s->associated, t, data);
    }
    return b;
}
