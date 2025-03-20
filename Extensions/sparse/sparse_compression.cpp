#include "sparse_compression.h"

static void
copy_sparsematrix_amatrix(psparsematrix src, uint roff, uint coff, pamatrix trg)
{
    assert(src->rows >= trg->rows + roff);
    assert(src->cols >= trg->cols + coff);

    clear_amatrix(trg);
    for (uint i = 0; i < trg->rows; i++)
    {
        uint ii = i + roff;
        for (uint k = src->row[ii]; k < src->row[ii + 1]; k++)
        {
            uint jj = src->col[k];

            if (jj >= coff + trg->cols)
            {
                // No more non-zero entries in this row
                // of the block, assuming ascending order
                // of the column indices
                break;
            }

            if (jj >= coff)
            {
                uint j = jj - coff;
                trg->a[i + j * trg->rows] = src->coeff[k];
            }
        }
    }
}

static void
copy_sparsematrix_rkmatrix(psparsematrix src, uint roff, uint coff, prkmatrix trg)
{
    uint rows, cols;
    rows = getrows_rkmatrix(trg);
    cols = getcols_rkmatrix(trg);

    assert(src->rows >= roff + rows);
    assert(src->cols >= coff + cols);
    // Quick exit, if matrix has no non-zero entries
    // Get number of non-zero rows of src
    // This is set as rank
    uint rank = 0;

    for (uint i = roff; i < roff + rows; i++)
    {
        uint start = src->row[i];
        uint end = src->row[i + 1];

        for (uint k = start; k < end; k++)
        {
            uint j = src->col[k];
            // Assume ascending order of the column indices
            if (j >= coff + cols)
            {
                break;
            }

            if (j >= coff && j < coff + cols)
            {
                rank++;
                break;
            }
        }
    }

    setrank_rkmatrix(trg, rank);

    pamatrix B = &trg->B;
    pamatrix A = &trg->A;
    clear_amatrix(A);
    clear_amatrix(B);

    // Copy non-zero rows of source matrix into B
    uint ell = 0;
    for (uint i = 0; i < rows; i++)
    {
        uint ii = i + roff;
        bool nzr = false;
        for (uint k = src->row[ii]; k < src->row[ii + 1]; k++)
        {
            uint jj = src->col[k];
            if (jj >= coff + cols)
            {
                // No non-zero entry left in this block-row,
                // assuming ascending order of the column indices
                break;
            }

            if (jj >= coff)
            {
                uint j = jj - coff;
                B->a[j + ell * cols] = src->coeff[k];
                nzr = true;
            }
        }

        if (nzr)
        {
            A->a[i + ell * rows] = 1.0;
            ell++;
        }
    }
}

static void
fastcopy_sparsematrix_to_hmatrix(psparsematrix src, uint roff, uint coff, phmatrix trg)
{
    assert(src->rows >= roff + trg->rc->size);
    assert(src->cols >= coff + trg->cc->size);
    // Matrix is subdevided
    if (trg->son != NULL)
    {
        // Assume that the source matrix is ordered like target matrix.
        // Column indices are assumed to be in ascending order.
        uint subcoff, subroff;
        subcoff = 0;
        subroff = 0;
        for (uint i = 0; i < trg->rsons; i++)
        {
            uint rows = trg->rc->son[i]->size;
            for (uint j = 0; j < trg->csons; j++)
            {
                uint cols = trg->cc->son[j]->size;

                fastcopy_sparsematrix_to_hmatrix(src, roff + subroff, coff + subcoff, trg->son[i + j * trg->rsons]);
                subcoff += cols;
            }
            subcoff = 0;
            subroff += rows;
        }
    }
    // dense block
    else if (trg->f != NULL)
    {
        copy_sparsematrix_amatrix(src, roff, coff, trg->f);
    }
    else
    {
        assert(trg->r != NULL);

        copy_sparsematrix_rkmatrix(src, roff, coff, trg->r);
    }
}

static inline psparsematrix
reorder(psparsematrix sp, uint *rowinv, uint *col)
{
    psparsematrix rsp;
    rsp = new_raw_sparsematrix(sp->rows, sp->cols, sp->nz);

    /* Reorder sparse matrix according to block structure of hm */
    uint ne = 0;
    for (uint i = 0; i < rsp->rows; i++)
    {
        rsp->row[i] = ne;
        /* Use identity for pinv in case of NULL pointer. */
        /* The pinv[i]-th row of A becomes the i-th row of B. */
        uint ii = rowinv[i];
        rsp->row[i] = ne;

        for (uint k = sp->row[ii]; k < sp->row[ii + 1]; k++)
        {
            rsp->coeff[ne] = sp->coeff[k];
            /* Use identity for q in case of NULL pointer. */
            /* The j-th column of A becomes the q[j]-th column of B. */
            rsp->col[ne] = col[sp->col[k]];
            /* Sort column indices. */
            for (uint h = ne; h > rsp->row[i]; h--)
            {
                if (rsp->col[h - 1] > rsp->col[h])
                {
                    uint j_tmp = rsp->col[h];
                    rsp->col[h] = rsp->col[h - 1];
                    rsp->col[h - 1] = j_tmp;

                    double coeff_tmp = rsp->coeff[h];
                    rsp->coeff[h] = rsp->coeff[h - 1];
                    rsp->coeff[h - 1] = coeff_tmp;
                }
            }
            ne += 1;
        }
    }
    assert(ne == sp->nz);
    rsp->row[rsp->rows] = ne;

    return rsp;
}

psparsematrix
reorder_sparsematrix(psparsematrix sp, uint *rowinv, uint *col)
{
    return reorder(sp, rowinv, col);
}

void copy_sparsematrix_to_hmatrix(psparsematrix sp, phmatrix hm)
{
    psparsematrix rsp;
    uint *rowinv = hm->rc->idx; // Inverse row permutation
    uint *colinv = hm->cc->idx; // Inverse column permutation
    uint *col;                  // Column permutation

    assert(hm->rc->size == sp->rows);
    assert(hm->cc->size == sp->cols);

    col = new uint[sp->cols];

    for (uint i = 0; i < sp->cols; i++)
    {
        col[colinv[i]] = i;
    }

    rsp = reorder(sp, rowinv, col);

    delete[] col;

    // Use reordered matrix to efficiently copy subblocks
    fastcopy_sparsematrix_to_hmatrix(rsp, 0, 0, hm);
    del_sparsematrix(rsp);
}