#include <string.h>
#include "oseencluster.h"

#ifdef USE_NETCDF
#include <stdio.h>
#include <netcdf.h>
#endif


inline void set_flags(uint *idx, uint size, psparsematrix sp, uint *flag)
{
  for(uint i = 0; i < size; i++)
  {
    uint ii = idx[i];
    for(uint k = sp->row[ii]; k < sp->row[ii+1]; k++)
    {
      flag[sp->col[k]] = 1;
    }
  }
}

inline void reset_flag(uint flag_size, uint *flag)
{
  for(uint i = 0; i < flag_size; i++)
  {
    flag[i] = 0;
  }
}

inline void swap_indices(uint *i, uint *j)
{
  uint tmp = *i;
  *i = *j;
  *j = tmp;
}

inline pcluster build_adaptive_interface_cluster(pclustergeometry cg, uint size,
                                                 uint *idx, uint clf, uint dim,
                                                 uint levelint)
{
  pcluster c;

  uint size0, size1;
  uint i, j, direction;
  real a, m;

  size0 = 0;
  size1 = 0;
  update_point_bbox_clustergeometry(cg, size, idx);

  if (size > clf)
  {
    /* compute the direction of partition */
    direction = 0;
    a = cg->hmax[0] - cg->hmin[0];
    for (j = 1; j < cg->dim; j++)
    {
      m = cg->hmax[j] - cg->hmin[j];
      if (a < m)
      {
        a = m;
        direction = j;
      }
    }

    if (levelint % dim)
    {
      levelint++;

      /* build sons */
      if (a > 0.0)
      {
        m = (cg->hmax[direction] + cg->hmin[direction]) / 2.0;
        size0 = 0;
        size1 = 0;

        for (i = 0; i < size; i++)
        {
          if (cg->x[idx[i]][direction] < m)
          {
            j = idx[i];
            idx[i] = idx[size0];
            idx[size0] = j;
            size0++;
          }
          else
          {
            size1++;
          }
        }
        if (size0 > 0)
        {
          if (size1 > 0)
          {
            c = new_cluster(size, idx, 2, cg->dim);

            c->son[0] = build_adaptive_interface_cluster(cg, size0, idx, clf,
                                                         dim, levelint);
            c->son[1] = build_adaptive_interface_cluster(cg, size1, idx + size0,
                                                         clf, dim, levelint);
            update_bbox_cluster(c);
          }
          else
          {
            assert(size0 > 0);
            assert(size1 == 0);
            c = new_cluster(size, idx, 1, cg->dim);
            c->son[0] = build_adaptive_interface_cluster(cg, size, idx, clf,
                                                         dim, levelint);
            update_bbox_cluster(c);
          }
        }
        else
        {
          assert(size0 == 0);
          assert(size1 > 0);
          c = new_cluster(size, idx, 1, cg->dim);
          c->son[0] = build_adaptive_interface_cluster(cg, size, idx, clf, dim,
                                                       levelint);
          update_bbox_cluster(c);
        }
      }
      else
      {
        assert(a == 0.0);
        c = new_cluster(size, idx, 0, cg->dim);
        update_support_bbox_cluster(cg, c);
      }
    }

    else
    {
      levelint++;
      c = new_cluster(size, idx, 1, cg->dim);
      c->son[0] =
          build_adaptive_interface_cluster(cg, size, idx, clf, dim, levelint);
      update_bbox_cluster(c);
    }
  }

  else
  {
    /* size <= clf */
    c = new_cluster(size, idx, 0, cg->dim);
    update_support_bbox_cluster(cg, c);
  }

  c->type = 2;
  update_cluster(c);

  return c;
}

pcluster build_simsub_interface_cluster(pclustergeometry cg, uint size,
                                        uint *idx, uint clf, uint direction,
                                        uint levelint)
{
  pcluster c;

  if (size > clf)
  {
    if (levelint % cg->dim)
    {
      levelint++;

      real a = cg->hmax[direction] - cg->hmin[direction];
      if (a > 0.0)
      {
        real m = (cg->hmax[direction] + cg->hmin[direction]) * 0.5;
        uint size0, size1;
        size0 = 0;
        size1 = 0;

        // Split index set
        for (uint i = 0; i < size; i++)
        {
          if (cg->x[idx[i]][direction] < m)
          {
            uint j = idx[i];
            idx[i] = idx[size0];
            idx[size0] = j;
            size0++;
          }
          else
            size1++;
        }

        if (size0 > 0)
        {
          if (size1 > 0)
          {
            c = new_cluster(size, idx, 2, cg->dim);
            real a, b;
            a = cg->hmin[direction];
            b = cg->hmax[direction];

            cg->hmax[direction] = m;
            c->son[0] = build_simsub_interface_cluster(cg, size0, idx, clf, (direction + 1) % cg->dim, levelint);
            cg->hmax[direction] = b;

            cg->hmin[direction] = m;
            c->son[1] = build_simsub_interface_cluster(cg, size1, idx + size0, clf, (direction + 1) % cg->dim, levelint);
            cg->hmin[direction] = a;

            update_bbox_cluster(c);
          }
          else
          {
            c = new_cluster(size, idx, 1, cg->dim);
            real b = cg->hmax[direction];

            cg->hmax[direction] = m;
            c->son[0] = build_simsub_interface_cluster(cg, size, idx, clf, (direction + 1) % cg->dim, levelint);
            cg->hmin[direction] = b;

            update_bbox_cluster(c);
          }
        }
        else
        {
          assert(size1 > 0);
          c = new_cluster(size, idx, 1, cg->dim);

          real a = cg->hmin[direction];

          cg->hmin[direction] = m;
          c->son[0] = build_simsub_interface_cluster(cg, size, idx, clf, (direction + 1) % cg->dim, levelint);
          cg->hmin[direction] = a;

          update_bbox_cluster(c);
        }
      }
      else
      {
        assert(a == 0.0);
        c = new_cluster(size, idx, 0, cg->dim);
        update_support_bbox_cluster(cg, c);
      }
    }
    else
    {
      levelint++;
      c = new_cluster(size, idx, 1, cg->dim);
      c->son[0] = build_simsub_interface_cluster(cg, size, idx, clf, (direction + 1) % cg->dim, levelint);
      update_bbox_cluster(c);
    }
  }
  else
  {
    c = new_cluster(size, idx, 0, cg->dim);
    update_support_bbox_cluster(cg, c);
  }

  c->type = 2;
  update_cluster(c);

  return c;
}

pcluster build_simsub_dd_cluster(pclustergeometry cg, uint size,
                                 uint *idx, uint clf, psparsematrix sp,
                                 uint direction, uint *flag)
{
  pcluster c;

  if (size > clf)
  {
    real a = cg->hmax[direction] - cg->hmin[direction];
    if (a > 0.0)
    {
      real m = (cg->hmax[direction] + cg->hmin[direction]) * 0.5;
      uint size0, size1, size2;

      size0 = 0;
      size1 = 0;
      size2 = 0;

      // Split index set into two parts
      for (uint i = 0; i < size; i++)
      {
        if (cg->x[idx[i]][direction] < m)
        {
          uint j = idx[i];
          idx[i] = idx[size0];
          idx[size0] = j;
          size0++;
        }
        else
          size1++;
      }

      for (uint i = 0; i < size0; i++)
        flag[idx[i]] = 1;

      // Form Interface cluster
      uint i = size0;
      while (i < size0 + size1)
      {
        bool inter = false;
        for (uint j = sp->row[idx[i]]; j < sp->row[idx[i] + 1]; j++)
        {
          if (flag[sp->col[j]] == 1)
          {
            uint tmp = idx[i];
            idx[i] = idx[size - size2 - 1];
            idx[size - size2 - 1] = tmp;
            size2++;
            size1--;
            inter = true;
            break;
          }
        }
        if (inter == false)
          i++;
      }

      // Reset flag array
      for (uint i = 0; i < size0; i++)
        flag[idx[i]] = 0;

      // Build sons
      if (size0 > 0)
        if (size1 > 0)
          if (size2 > 0)
          {
            c = new_cluster(size, idx, 3, cg->dim);

            real a, b;
            a = cg->hmin[direction];
            b = cg->hmax[direction];

            cg->hmax[direction] = m;
            c->son[0] = build_simsub_dd_cluster(cg, size0, idx, clf, sp, (direction + 1) % cg->dim, flag);
            cg->hmin[direction] = m;
            c->son[2] = build_simsub_interface_cluster(cg, size2, idx + size0 + size1,
                                                       clf, (direction + 1) % cg->dim, 1);
            cg->hmax[direction] = b;
            c->son[1] = build_simsub_dd_cluster(cg, size1, idx + size0, clf, sp, (direction + 1) % cg->dim, flag);
            cg->hmin[direction] = a;

            update_bbox_cluster(c);
          }
          else
          {
            c = new_cluster(size, idx, 2, cg->dim);
            real a, b;
            a = cg->hmin[direction];
            b = cg->hmax[direction];

            cg->hmax[direction] = m;
            c->son[0] = build_simsub_dd_cluster(cg, size0, idx, clf, sp, (direction + 1) % cg->dim, flag);
            cg->hmax[direction] = b;

            cg->hmin[direction] = m;
            c->son[1] = build_simsub_dd_cluster(cg, size1, idx + size0, clf, sp, (direction + 1) % cg->dim, flag);
            cg->hmin[direction] = a;

            update_bbox_cluster(c);
          }
        else
        {
          c = new_cluster(size, idx, 1, cg->dim);
          c->son[0] = build_simsub_dd_cluster(cg, size, idx, clf, sp, (direction + 1) % cg->dim, flag);
          update_bbox_cluster(c);
        }
      else
      {
        c = new_cluster(size, idx, 1, cg->dim);
        c->son[0] = build_simsub_dd_cluster(cg, size, idx, clf, sp, (direction + 1) % cg->dim, flag);
        update_bbox_cluster(c);
      }
    }
    else
    {
      c = new_cluster(size, idx, 1, cg->dim);
      c->son[0] = build_simsub_dd_cluster(cg, size, idx, clf, sp, (direction + 1) % cg->dim, flag);
      update_bbox_cluster(c);
    }
  }
  else
  {
    c = new_cluster(size, idx, 0, cg->dim);
    update_support_bbox_cluster(cg, c);
  }

  c->type = 1;
  update_cluster(c);

  return c;
}

pcluster build_simsub_cluster_level(pclustergeometry cg, uint size,
                                    uint *idx, uint clf, uint direction)
{
  pcluster c;

  if (clf < size)
  {
    real width = cg->hmax[direction] - cg->hmin[direction];
    if (width > 0.0)
    {
      real mean = (cg->hmax[direction] + cg->hmin[direction]) * 0.5;
      uint size0, size1;

      size0 = 0;
      size1 = 0;
      // Devide index set according to geometrical information
      for (uint i = 0; i < size; i++)
      {
        if (cg->x[idx[i]][direction] < mean)
        {
          // Swap entry i to the beginning of idx
          uint tmp = idx[i];
          idx[i] = idx[size0];
          idx[size0] = tmp;
          // Increase number of nodes in the first half
          size0++;
        }
        else
        {
          // Increase number of nodes in the second half
          size1++;
        }
      }

      assert(size == size0 + size1);

      // Build sons
      if (size0 > 0)
      {
        if (size1 > 0)
        {
          c = new_cluster(size, idx, 2, cg->dim);
          real a, b;

          a = cg->hmin[direction];
          b = cg->hmax[direction];

          cg->hmax[direction] = mean;
          c->son[0] = build_simsub_cluster_level(cg, size0, idx, clf, (direction + 1) % cg->dim);
          cg->hmax[direction] = b;

          cg->hmin[direction] = mean;
          c->son[1] = build_simsub_cluster_level(cg, size1, idx + size0, clf, (direction + 1) % cg->dim);
          cg->hmin[direction] = a;

          update_bbox_cluster(c);
        }
        else
        {
          c = new_cluster(size, idx, 1, cg->dim);

          real b = cg->hmax[direction];
          cg->hmax[direction] = mean;
          c->son[0] = build_simsub_cluster_level(cg, size, idx, clf, direction);
          cg->hmax[direction] = b;

          update_bbox_cluster(c);
        }
      }
      else
      {
        c = new_cluster(size, idx, 1, cg->dim);

        real a = cg->hmin[direction];
        cg->hmin[direction] = mean;
        c->son[0] = build_simsub_cluster_level(cg, size, idx, clf, (direction + 1) % cg->dim);
        cg->hmin[direction] = a;

        update_bbox_cluster(c);
      }
    }
    else
    {
      assert(width == 0.0);
      c = new_cluster(size, idx, 1, cg->dim);
      c->son[0] = build_simsub_cluster_level(cg, size, idx, clf, (direction + 1) % cg->dim);
      update_support_bbox_cluster(cg, c);
    }
  }
  else
  {
    c = new_cluster(size, idx, 0, cg->dim);

    update_support_bbox_cluster(cg, c);
  }

  update_cluster(c);
  return c;
}

uint *sort_indices_flag(uint size, uint *idx, uint *flag)
{
  // Get the maximal flag (either 1 or 2)
  uint max = 0;
  for (uint i = 0; i < size; i++)
  {
    uint ii = idx[i];
    max = (max < flag[ii] ? flag[ii] : max);
  }

  uint size0 = 0;
  uint size1 = 0;
  uint size2 = 0;
  if (max >= 1)
  {
    // Sort index set
    for (uint i = 0; i < size; i++)
    {
      if (flag[idx[i]] == 0)
      {
        uint tmp = idx[i];
        idx[i] = idx[size0];
        idx[size0] = tmp;
        size0++;
      }
      else
      {
        size1++;
      }
    }
  }
  if (max == 2)
  {
    uint j = size0;
    while (j < size - size2)
    {
      if (flag[idx[j]] == 2)
      {
        uint tmp = idx[j];
        idx[j] = idx[size - size2 - 1];
        idx[size - size2 - 1] = tmp;
        size2++;
        size1--;
      }
      else
        j++;
    }
  }

  uint *sizes = new uint[3];
  sizes[0] = size0;
  sizes[1] = size1;
  sizes[2] = size2;

  return sizes;
}

pcluster build_coupled_dd_cluster(pcluster cluster_p, pclustergeometry cg, uint size,
                                  uint *idx, uint clf, psparsematrix graph, uint *flag)
{
  pcluster c;

  if (cluster_p->sons > 1)
  {
    pcluster son;
    assert(cluster_p->sons == 2);

    // Set flag for nodes connected to the first cluster
    son = cluster_p->son[1];
    for (uint i = 0; i < son->size; i++)
    {
      uint ii = son->idx[i];
      for (uint k = graph->row[ii]; k < graph->row[ii + 1]; k++)
        flag[graph->col[k]] = 1;
    }

    uint size0 = 0;
    uint size1 = 0;
    uint size2 = 0;
    for (uint i = 0; i < size; i++)
    {
      if (flag[idx[i]] == 0)
      {
        uint tmp = idx[i];
        idx[i] = idx[size0];
        idx[size0] = tmp;
        size0++;
      }
      else
        size1++;
    }

    // Reset flag
    for (uint i = 0; i < graph->cols; i++)
      flag[i] = 0;

    son = cluster_p->son[0];
    for (uint i = 0; i < son->size; i++)
    {
      uint ii = son->idx[i];
      for (uint k = graph->row[ii]; k < graph->row[ii + 1]; k++)
        flag[graph->col[k]] = 1;
    }

    uint j = size0;
    while (j < size0 + size1)
    {
      if (flag[idx[j]] == 1)
      {
        uint tmp = idx[j];
        idx[j] = idx[size - size2 - 1];
        idx[size - size2 - 1] = tmp;
        size1--;
        size2++;
      }
      else
      {
        j++;
      }
    }

    // Reset flag
    for (uint i = 0; i < graph->cols; i++)
      flag[i] = 0;

    if (size0 > 0)
    {
      if (size1 > 0)
      {
        if (size2 > 0)
        {
          c = new_cluster(size, idx, 3, cg->dim);
          c->son[0] = build_coupled_dd_cluster(cluster_p->son[0], cg, size0, idx, clf, graph, flag);
          c->son[1] = build_coupled_dd_cluster(cluster_p->son[1], cg, size1, idx + size0, clf, graph, flag);
          c->son[2] = build_adaptive_interface_cluster(cg, size2, idx + size0 + size1, clf, cg->dim, 1);

          update_bbox_cluster(c);
        }
        else
        {
          c = new_cluster(size, idx, 2, cg->dim);
          c->son[0] = build_coupled_dd_cluster(cluster_p->son[0], cg, size0, idx, clf, graph, flag);
          c->son[1] = build_coupled_dd_cluster(cluster_p->son[1], cg, size1, idx + size0, clf, graph, flag);

          update_bbox_cluster(c);
        }
      }
      else
      {
        if (size2 > 0)
        {
          c = new_cluster(size, idx, 2, cg->dim);
          c->son[0] = build_coupled_dd_cluster(cluster_p->son[0], cg, size0, idx, clf, graph, flag);
          c->son[1] = build_adaptive_interface_cluster(cg, size2, idx + size0, clf, cg->dim, 1);

          update_bbox_cluster(c);
        }
        else
        {

          c = new_cluster(size, idx, 1, cg->dim);
          c->son[0] = build_coupled_dd_cluster(cluster_p->son[0], cg, size, idx, clf, graph, flag);

          update_bbox_cluster(c);
        }
      }
    }
    else
    {
      if (size1 > 0)
      {
        if (size2 > 0)
        {
          c = new_cluster(size, idx, 2, cg->dim);
          c->son[0] = build_coupled_dd_cluster(cluster_p->son[0], cg, size1, idx, clf, graph, flag);
          c->son[1] = build_adaptive_interface_cluster(cg, size2, idx + size1, clf, cg->dim, 1);

          update_bbox_cluster(c);
        }
        else
        {
          c = new_cluster(size, idx, 1, cg->dim);
          c->son[0] = build_coupled_dd_cluster(cluster_p->son[0], cg, size, idx, clf, graph, flag);

          update_bbox_cluster(c);
        }
      }
      else
      {
        c = new_cluster(size, idx, 1, cg->dim);
        c->son[0] = build_adaptive_interface_cluster(cg, size, idx, clf, cg->dim, 1);

        update_bbox_cluster(c);
      }
    }
  }
  else if (cluster_p->sons == 1)
  {
    // cluster_p has only one son
    c = new_cluster(size, idx, 1, cg->dim);

    c->son[0] = build_coupled_dd_cluster(cluster_p->son[0], cg, size, idx, clf, graph, flag);

    update_bbox_cluster(c);
  }
  else
  {
    c = new_cluster(size, idx, 0, cg->dim);

    update_support_bbox_cluster(cg, c);
  }

  update_cluster(c);

  c->type = 1;
  c->associated = cluster_p;
  cluster_p->associated = c;

  return c;
}

pcluster build_coupled_dd_cluster_v2(pcluster cluster_p, pclustergeometry cg, uint size, uint *idx, 
                                     uint clf, psparsematrix sp_coupling, psparsematrix sp_vel, uint *flag)
{
  pcluster c;

  if (cluster_p->sons > 1)
  {
    pcluster son;
    assert(cluster_p->sons == 2);

    // Set flag for nodes connected to the first cluster
    son = cluster_p->son[1];
    for (uint i = 0; i < son->size; i++)
    {
      uint ii = son->idx[i];
      for (uint k = sp_coupling->row[ii]; k < sp_coupling->row[ii + 1]; k++)
        flag[sp_coupling->col[k]] = 1;
    }

    uint size0 = 0;
    uint size1 = 0;
    for (uint i = 0; i < size; i++)
    {
      if (flag[idx[i]] == 0)
      {
        uint tmp = idx[i];
        idx[i] = idx[size0];
        idx[size0] = tmp;
        size0++;
      }
      else
        size1++;
    }

    // Reset flag
    for (uint i = 0; i < sp_coupling->cols; i++)
      flag[i] = 0;

    son = cluster_p->son[0];
    for (uint i = 0; i < son->size; i++)
    {
      uint ii = son->idx[i];
      for (uint k = sp_coupling->row[ii]; k < sp_coupling->row[ii + 1]; k++)
        flag[sp_coupling->col[k]] = 1;
    }

    uint size2 = 0;
    uint j = size0;
    while (j < size0 + size1)
    {
      if (flag[idx[j]] == 1)
      {
        uint tmp = idx[j];
        idx[j] = idx[size - size2 - 1];
        idx[size - size2 - 1] = tmp;
        size1--;
        size2++;
      }
      else
      {
        j++;
      }
    }

    // Reset flag
    for (uint i = 0; i < sp_coupling->cols; i++)
      flag[i] = 0;

    // Split interface into three parts
    // One connected to the left domain cluster
    // One connected to the right domain cluster
    // One connected to neither of them 
    for(uint i = 0; i < size0; i++)
    {
      uint ii = idx[i];
      for(uint k = sp_vel->row[ii]; k < sp_vel->row[ii+1]; k++)
      {
        flag[sp_vel->col[k]] = 1;
      }
    }
    
    uint current_size = size2;

    uint size3 = size2;
    size2 = 0;
    uint *current = idx + size0 + size1;

    for(uint i = 0; i < current_size; i++)
    {
      if(flag[current[i]] == 1)
      {
        uint tmp = current[i];
        current[i] = current[size2];
        current[size2] = tmp;
        size2++;
        size3--;
      }
    }

    // Reset flag
    for (uint i = 0; i < sp_coupling->cols; i++)
      flag[i] = 0;

    for(uint i = size0; i < size0 + size1; i++)
    {
      uint ii = idx[i];
      for(uint k = sp_vel->row[ii]; k < sp_vel->row[ii+1]; k++)
      {
        flag[sp_vel->col[k]] = 1;
      }
    }

    current_size = size3;
    current = current + size2;
    uint size4 = size3;
    size3 = 0;
    for(uint i = 0; i < current_size; i++)
    {
      if(flag[current[i]] == 0)
      {
        uint tmp = current[i];
        current[i] = current[size3];
        current[size3] = tmp;
        size3++;
        size4--;
      }
    }

    // Reset flag
    for (uint i = 0; i < sp_coupling->cols; i++)
      flag[i] = 0;

    assert(size0 + size1 + size2 + size3 + size4 == size);

    uint *sizes;
    sizes = new uint[5];
    sizes[0] = size0;
    sizes[1] = size1;
    sizes[2] = size2;
    sizes[3] = size3;
    sizes[4] = size4;

    // Initialize array for sons with NULL's
    pcluster *sons;
    sons = new pcluster[5];
    for(uint i = 0; i < 5; i++)
    {
      sons[i] = NULL;
    }
    uint son_count = 0;

    // Build domain cluster sons
    uint * sidx = idx;
    for(uint i = 0; i < 2; i++)
    {
      if(sizes[i] > 0)
      {
        sons[i] = build_coupled_dd_cluster_v2(cluster_p->son[i], cg, sizes[i], sidx, clf, sp_coupling, sp_vel, flag);
        sidx += sizes[i];

        son_count++;
      }
    }

    // Build left interface son 
    if(sizes[2] > 0)
    {
      sons[2] = build_adaptive_interface_cluster(cg, sizes[2], sidx, clf, cg->dim, 1);
      sidx += sizes[2];

      sons[2]->associated = sons[0];
      sons[2]->type = 3;


      son_count++;
    }

    // Build middle interface son
    if(sizes[3] > 0)
    {
      sons[3] = build_adaptive_interface_cluster(cg, sizes[3], sidx, clf, cg->dim, 1);
      sidx += sizes[3];

      sons[3]->type = 4;

      son_count++;
    }

    // Build right interface son
    if(sizes[4] > 0)
    {
      sons[4] = build_adaptive_interface_cluster(cg, sizes[4], sidx, clf, cg->dim, 1);
      sidx += sizes[4];

      sons[4]->associated = sons[1];
      sons[4]->type = 3;

      son_count++;
    }

    c = new_cluster(size, idx, son_count, cg->dim);
    uint k = 0;
    for(uint i = 0; i < 5; i++)
    {
      if(sons[i] != NULL)
      {
        c->son[k] = sons[i];
        k++;
      }
    }
    assert(k == son_count);

    update_bbox_cluster(c);

    delete[] sizes;
    delete[] sons;

  }
  else if (cluster_p->sons == 1)
  {
    // cluster_p has only one son
    c = new_cluster(size, idx, 1, cg->dim);

    c->son[0] = build_coupled_dd_cluster(cluster_p->son[0], cg, size, idx, clf, sp_coupling, flag);

    update_bbox_cluster(c);
  }
  else
  {
    c = new_cluster(size, idx, 0, cg->dim);

    update_support_bbox_cluster(cg, c);
  }

  update_cluster(c);

  c->type = 1;
  c->associated = cluster_p;
  cluster_p->associated = c;

  return c;
}


#ifdef USE_NETCDF
static void
write_count(pccluster t, size_t *clusters, size_t *coeffs)
{
  uint i;

  /* Increase cluster counter */
  (*clusters)++;

  /* Add number of coefficients of transfer matrix */
  (*coeffs) += 2 * t->dim;

  /* Handle sons */
  for (i = 0; i < t->sons; i++)
    write_count(t->son[i], clusters, coeffs);
}

static void
write_cdf(pccluster t,
          size_t clusters, size_t coeffs,
          size_t *clusteridx, size_t *coeffidx,
          int nc_file, int nc_sons, int nc_size,
          int nc_type, int nc_coeff)
{
  size_t start, count;
  ptrdiff_t stride;
  int val, result;
  uint i;

  assert(*clusteridx <= clusters);

  /* Write number of sons to nc_sons[*clusteridx] */
  start = *clusteridx;
  count = 1;
  stride = 1;
  val = t->sons;
  result = nc_put_vars(nc_file, nc_sons, &start, &count, &stride, &val);
  assert(result == NC_NOERR);

  /* Write size of cluster to nc_size[*clusteridx] */
  val = t->size;
  result = nc_put_vars(nc_file, nc_size, &start, &count, &stride, &val);
  assert(result == NC_NOERR);

  val = t->type;
  result = nc_put_vars(nc_file, nc_type, &start, &count, &stride, &val);
  assert(result == NC_NOERR);

  /* Increase cluster index */
  (*clusteridx)++;

  /* Handle sons */
  for (i = 0; i < t->sons; i++)
    write_cdf(t->son[i], clusters, coeffs, clusteridx, coeffidx,
              nc_file, nc_sons, nc_size, nc_type, nc_coeff);

  /* Write bounding box */
  start = *coeffidx;
  assert(start + 2 * t->dim <= coeffs);
  count = t->dim;
  result = nc_put_vars(nc_file, nc_coeff, &start, &count, &stride, t->bmin);
  assert(result == NC_NOERR);
  start += t->dim;

  result = nc_put_vars(nc_file, nc_coeff, &start, &count, &stride, t->bmax);
  assert(result == NC_NOERR);
  start += t->dim;
  (*coeffidx) = start;
}

void write_cdf_cluster_grid(pccluster t, pclustergeometry cg, const char *name)
{
  size_t clusters, clusteridx;
  size_t coeffs, coeffidx;
  int nc_file, nc_sons, nc_size, nc_idx, nc_coeff;
  int nc_clusters, nc_coeffs, nc_totalsize, nc_dim;
  int nc_x, nc_ndofs, nc_totalcoords;
  int nc_type;
  int result;

  /* Count number of clusters and coefficients */
  clusters = 0;
  coeffs = 0;
  write_count(t, &clusters, &coeffs);

  /* Create NetCDF file */
  result = nc_create(name, NC_64BIT_OFFSET, &nc_file);
  assert(result == NC_NOERR);

  /* Define "clusters" dimension */
  result = nc_def_dim(nc_file, "clusters", clusters, &nc_clusters);
  assert(result == NC_NOERR);

  /* Define "coeffs" dimension */
  result = nc_def_dim(nc_file, "coeffs", coeffs, &nc_coeffs);
  assert(result == NC_NOERR);

  /* Define "totalsize" dimension */
  result = nc_def_dim(nc_file, "totalsize", t->size, &nc_totalsize);
  assert(result == NC_NOERR);

  /* Define "dim" dimension */
  result = nc_def_dim(nc_file, "dim", t->dim, &nc_dim);
  assert(result == NC_NOERR);

  result = nc_def_dim(nc_file, "ndofs", cg->nidx, &nc_ndofs);
  assert(result == NC_NOERR);

  result = nc_def_dim(nc_file, "totalcoords", cg->nidx * 3, &nc_totalcoords);
  assert(result == NC_NOERR);

  /* Define "sons" variable */
  result = nc_def_var(nc_file, "sons", NC_INT, 1, &nc_clusters, &nc_sons);
  assert(result == NC_NOERR);

  /* Define "size" variable */
  result = nc_def_var(nc_file, "size", NC_INT, 1, &nc_clusters, &nc_size);
  assert(result == NC_NOERR);

  /* Define "idx" variable */
  result = nc_def_var(nc_file, "idx", NC_INT, 1, &nc_totalsize, &nc_idx);
  assert(result == NC_NOERR);

  /* Define "coeff" variable */
  result = nc_def_var(nc_file, "coeff", NC_DOUBLE, 1, &nc_coeffs, &nc_coeff);
  assert(result == NC_NOERR);

  /* Define "gridpoints" variable for grid points */
  result = nc_def_var(nc_file, "gridpoints", NC_DOUBLE, 1, &nc_totalcoords, &nc_x);
  assert(result == NC_NOERR);

  result = nc_def_var(nc_file, "types", NC_INT, 1, &nc_clusters, &nc_type);
  assert(result == NC_NOERR);

  /* Finish NetCDF define mode */
  result = nc_enddef(nc_file);
  assert(result == NC_NOERR);

  /* Write index to NetCDF variable */
  result = nc_put_var(nc_file, nc_idx, t->idx);

  /* Write grid points to NetCDF variable */
  real *xvec = new real[cg->nidx * 3];
  for (uint i = 0; i < cg->nidx; i++)
    for (uint j = 0; j < 3; j++)
      xvec[j + i * 3] = cg->x[i][j];

  result = nc_put_var(nc_file, nc_x, xvec);

  /* Write coefficiens to NetCDF variables */
  clusteridx = 0;
  coeffidx = 0;
  write_cdf(t, clusters, coeffs, &clusteridx, &coeffidx,
            nc_file, nc_sons, nc_size, nc_type, nc_coeff);
  assert(clusteridx == clusters);
  assert(coeffidx == coeffs);

  /* Close file */
  result = nc_close(nc_file);
  assert(result == NC_NOERR);

  delete[] xvec;
}

static void
prefix_name(char *buf, int bufsize, const char *prefix, const char *name)
{
  if (prefix)
    snprintf(buf, bufsize, "%s_%s", prefix, name);
  else
    strncpy(buf, name, bufsize);
}
#endif
