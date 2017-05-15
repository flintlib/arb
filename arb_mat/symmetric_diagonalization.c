/*=============================================================================

    This file is part of ARB.

    ARB is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    ARB is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with ARB; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2016 Jonathan Bober, etc.

******************************************************************************/

#include "arb_mat.h"

struct _sortme
{
    arb_ptr p;
    int idx;
};

static int
_arb_cmp_for_sort(const void *a, const void *b)
{
    const struct _sortme *x = a;
    const struct _sortme *y = b;
    return arf_cmp(arb_midref(x->p), arb_midref(y->p));
}

/* midpoints of diagonal entries will be non-decreasing */
static void
_sort_decomposition(arb_mat_t D, arb_mat_t P)
{
    slong j, n;
    struct _sortme *s;

    n = arb_mat_nrows(P);
    s = flint_malloc(n * sizeof(struct _sortme));

    for (j = 0; j < n; j++)
    {
        s[j].p = arb_mat_entry(D, j, 0);
        s[j].idx = j;
    }

    qsort(s, n, sizeof(struct _sortme), _arb_cmp_for_sort);

    {
        slong i;
        arb_mat_t Pt, Dt;
        arb_mat_init(Pt, n, n);
        arb_mat_init(Dt, n, 1);
        arb_mat_set(Pt, P);
        arb_mat_set(Dt, D);
        for (j = 0; j < n; j++)
        {
            arb_set(arb_mat_entry(D, j, 0),
                    arb_mat_entry(Dt, s[j].idx, 0));
            for (i = 0; i < n; i++)
            {
                arb_set(arb_mat_entry(P, i, j),
                        arb_mat_entry(Pt, i, s[j].idx));
            }
        }
        arb_mat_clear(Pt);
        arb_mat_clear(Dt);
    }

    flint_free(s);
}

/* multiplies each column by the sign of its leading nonzero midpoint */
static void
_standardize_column_signs(arb_mat_t P)
{
    slong i, j;
    for (j = 0; j < arb_mat_ncols(P); j++)
    {
        int negate_col = 0;
        for (i = 0; i < arb_mat_nrows(P); i++)
        {
            int sgn = arf_sgn(arb_midref(arb_mat_entry(P, i, j)));
            if (sgn < 0)
            {
                negate_col = 1;
                break;
            }
            else if (sgn > 0)
            {
                break;
            }
        }
        if (negate_col)
        {
            for (i = 0; i < arb_mat_nrows(P); i++)
            {
                arb_neg(arb_mat_entry(P, i, j), arb_mat_entry(P, i, j));
            }
        }
    }
}

/* divides each column by its norm */
static void
_standardize_column_norms(arb_mat_t P, slong prec)
{
    slong i, j;
    arb_t t;
    arb_init(t);
    for (j = 0; j < arb_mat_ncols(P); j++)
    {
        arb_zero(t);
        for (i = 0; i < arb_mat_nrows(P); i++)
        {
            arb_srcptr z = arb_mat_entry(P, i, j);
            arb_addmul(t, z, z, prec);
        }
        arb_sqrtpos(t, t, prec);
        for (i = 0; i < arb_mat_nrows(P); i++)
        {
            arb_ptr z = arb_mat_entry(P, i, j);
            arb_div(z, z, t, prec);
        }
    }
    arb_clear(t);
}


static void
_arf_hypot(arf_t z, const arf_t x, const arf_t y, slong prec)
{
    if (arf_is_zero(y))
    {
        arf_abs(z, x);
    }
    else if (arf_is_zero(x))
    {
        arf_abs(z, y);
    }
    else
    {
        arf_t t;
        arf_init(t);
        arf_mul(t, x, x, prec + 4, ARB_RND);
        arf_mul(z, y, y, prec + 4, ARB_RND);
        arf_add(t, t, z, prec + 4, ARB_RND);
        arf_sqrt(z, t, prec, ARB_RND);
        arf_clear(t);
    }
}

static void
_arf_twobytwo_diag(arf_t u1, arf_t u2,
        const arf_t a, const arf_t b, const arf_t d, slong prec) {
    /*
    // Compute the orthogonal matrix that diagonalizes
    //
    //    A = [a b]
    //        [b d]
    //
    // This matrix will have the form
    //
    //    U = [cos x , -sin x]
    //        [sin x, cos x]
    //
    // where the diagonal matrix is U^t A U.
    // We set u1 = cos x, u2 = -sin x.
    */
    arf_t x;
    arf_t r;

    if (arf_is_zero(b))
    {
        arf_set_ui(u1, 1);
        arf_set_ui(u2, 0);
        return;
    }
    arf_init(x);
    arf_init(r);

    /* r = (a-d)/2 */
    arf_sub(r, a, d, prec, ARF_RND_NEAR);
    arf_mul_2exp_si(r, r, -1);

    /* u1 = sqrt(b^2 + ((a-d)/2)^2) */
    _arf_hypot(u1, b, r, prec);

    /* u1 = (d - a)/2 + sqrt(b^2 + ( (a-d)/2 )^2) */
    arf_sub(u1, u1, r, prec, ARF_RND_NEAR);

    /* x = sqrt(u1^2 + b^2) */
    _arf_hypot(x, u1, b, prec);
    arf_div(u2, u1, x, prec, ARF_RND_NEAR);
    arf_div(u1, b, x, prec, ARF_RND_NEAR);
    arf_neg(u1, u1);

    arf_clear(x);
    arf_clear(r);
}


int
_arb_mat_jacobi_diagonalization(arb_mat_t D, arb_mat_t P, const arb_mat_t A, slong prec) {

#define B(i,j) arb_mat_entry(B, i, j)
#define D(i) arb_mat_entry(D, i, 0)
#define P(i,j) arb_mat_entry(P, i, j)

    int unique_eigenvalues = 1;
    arf_struct *B1, *B2, *row_max;
    int *row_max_indices;
    int dim, i, j, k, l;
    int iter, finished;
    arb_mat_t B;
    arf_t x1;
    arf_t Gii, Gij, Gji, Gjj;

    if (!arb_mat_is_square(A))
    {
        flint_printf("Exception (_arb_mat_jacobi_diagonalization). "
                     "Non-square matrix.\n");
        abort();
    }

    dim = arb_mat_nrows(A);

    if (arb_mat_nrows(P) != dim || arb_mat_ncols(P) != dim ||
        arb_mat_nrows(D) != dim || arb_mat_ncols(D) != 1)
    {
        flint_printf("Exception (_arb_mat_jacobi_diagonalization). "
                     "Incompatible matrix dimensions.\n");
        abort();
    }

    if (dim == 0)
        return 0;

    if (dim == 1)
    {
        arb_mat_set(D, A);
        arb_mat_one(P);
        return 0;
    }
    arb_mat_init(B, dim, dim);

    B1 = flint_malloc(dim * sizeof(arf_struct));
    B2 = flint_malloc(dim * sizeof(arf_struct));
    row_max = flint_malloc((dim - 1) * sizeof(arf_struct));
    row_max_indices = flint_malloc((dim - 1) * sizeof(int));

    for (k = 0; k < dim; k++)
    {
        arf_init(B1+k);
        arf_init(B2+k);
    }
    for (k = 0; k < dim - 1; k++)
    {
        arf_init(row_max+k);
    }

    arf_init(x1);

    arf_init(Gii);
    arf_init(Gij);
    arf_init(Gji);
    arf_init(Gjj);

    arb_mat_set(B, A);
    arb_mat_one(P);

    for (i = 0; i < dim - 1; i++)
    {
        for (j = i + 1; j < dim; j++)
        {
            arf_abs(x1, arb_midref(B(i,j)));
            if (arf_cmp(row_max+i, x1) <= 0)
            {
                arf_set(row_max+i, x1);
                row_max_indices[i] = j;
            }
        }
    }

    finished = 0;
    iter = 0;
    while (!finished)
    {
        slong bound;

        iter++;

        arf_zero(x1);
        i = 0;
        j = 0;
        for (k = 0; k < dim - 1; k++)
        {
            if (arf_cmp(x1, row_max+k) < 0)
            {
                arf_set(x1, row_max+k);
                i = k;
            }
        }
        j = row_max_indices[i];

        bound = arf_abs_bound_lt_2exp_si(x1);

        if (bound < -prec * .9)
        {
            finished = 1;
            break;
        }

        _arf_twobytwo_diag(
                Gii, Gij,
                arb_midref(B(i,i)), arb_midref(B(i,j)), arb_midref(B(j,j)),
                2*prec);
        arf_neg(Gji, Gij);
        arf_set(Gjj, Gii);

        if (arf_is_zero(Gij))
        {                       /* If this happens, we're */
            finished = 1;       /* not going to do any better */
            break;              /* without increasing the precision. */
        }

        for (k = 0; k < dim; k++)
        {
            arf_mul(B1+k, Gii, arb_midref(B(i,k)), prec, ARF_RND_NEAR);
            arf_addmul(B1+k, Gji, arb_midref(B(j,k)), prec, ARF_RND_NEAR);

            arf_mul(B2+k, Gij, arb_midref(B(i,k)), prec, ARF_RND_NEAR);
            arf_addmul(B2+k, Gjj, arb_midref(B(j,k)), prec, ARF_RND_NEAR);
        }
        for (k = 0; k < dim; k++)
        {
            arf_set(arb_midref(B(i,k)), B1+k);
            arf_set(arb_midref(B(j,k)), B2+k);
        }

        for (k = 0; k < dim; k++)
        {
            arf_mul(B1+k, Gii, arb_midref(B(k,i)), prec, ARF_RND_NEAR);
            arf_addmul(B1+k, Gji, arb_midref(B(k,j)), prec, ARF_RND_NEAR);

            arf_mul(B2+k, Gij, arb_midref(B(k,i)), prec, ARF_RND_NEAR);
            arf_addmul(B2+k, Gjj, arb_midref(B(k,j)), prec, ARF_RND_NEAR);
        }
        for (k = 0; k < dim; k++)
        {
            arf_set(arb_midref(B(k,i)), B1+k);
            arf_set(arb_midref(B(k,j)), B2+k);
        }

        for (k = 0; k < dim; k++)
        {
            arf_mul(B1+k, Gii, arb_midref(P(k,i)), prec, ARF_RND_NEAR);
            arf_addmul(B1+k, Gji, arb_midref(P(k,j)), prec, ARF_RND_NEAR);

            arf_mul(B2+k, Gij, arb_midref(P(k,i)), prec, ARF_RND_NEAR);
            arf_addmul(B2+k, Gjj, arb_midref(P(k,j)), prec, ARF_RND_NEAR);
        }
        for (k = 0; k < dim; k++)
        {
            arf_set(arb_midref(P(k,i)), B1+k);
            arf_set(arb_midref(P(k,j)), B2+k);
        }

        /* Declare that the 2x2 diagonalization was successful. */
        arf_zero(arb_midref(B(i, j)));
        arf_zero(arb_midref(B(j, i)));

        if (i < dim - 1)
            arf_zero(row_max+i);
        if (j < dim - 1)
            arf_zero(row_max+j);

        /* Update the max in any row where the maximum
         * was in a column that changed. */
        for (k = 0; k < dim - 1; k++)
        {
            if (row_max_indices[k] == j || row_max_indices[k] == i)
            {
                arf_abs(row_max+k, arb_midref(B(k,k+1)));
                row_max_indices[k] = k+1;
                for (l = k+2; l < dim; l++)
                {
                    arf_abs(x1, arb_midref(B(k,l)));
                    if (arf_cmp(row_max+k, x1) < 0)
                    {
                        arf_set(row_max+k, x1);
                        row_max_indices[k] = l;
                    }
                }
            }
        }

        /* Update the max in the ith row. */
        for (k = i + 1; k < dim; k++)
        {
            arf_abs(x1, arb_midref(B(i, k)));
            if (arf_cmp(row_max+i, x1) < 0)
            {
                arf_set(row_max+i, x1);
                row_max_indices[i] = k;
            }
        }

        /* Update the max in the jth row. */
        for (k = j + 1; k < dim; k++)
        {
            arf_abs(x1, arb_midref(B(j, k)));
            if (arf_cmp(row_max+j, x1) < 0)
            {
                arf_set(row_max+j, x1);
                row_max_indices[j] = k;
            }
        }

        /* Go through column i to see if any of
         * the new entries are larger than the
         * max of their row. */
        for (k = 0; k < i; k++)
        {
            if (k == dim)
                continue;
            arf_abs(x1, arb_midref(B(k, i)));
            if (arf_cmp(row_max+k, x1) < 0)
            {
                arf_set(row_max+k, x1);
                row_max_indices[k] = i;
            }
        }

        /* And then column j. */
        for (k = 0; k < j; k++)
        {
            if (k == dim)
                continue;
            arf_abs(x1, arb_midref(B(k, j)));
            if (arf_cmp(row_max+k, x1) < 0)
            {
                arf_set(row_max+k, x1);
                row_max_indices[k] = j;
            }
        }
    }

    for (k = 0; k < dim; k++)
        arb_get_mid_arb(D(k), B(k, k));

    arf_clear(x1);
    arb_mat_clear(B);
    for (k = 0; k < dim; k++)
    {
        arf_clear(B1+k);
        arf_clear(B2+k);
    }
    for (k = 0; k < dim - 1; k++)
    {
        arf_clear(row_max+k);
    }
    arf_clear(Gii);
    arf_clear(Gij);
    arf_clear(Gji);
    arf_clear(Gjj);
    flint_free(B1);
    flint_free(B2);
    flint_free(row_max);
    flint_free(row_max_indices);

    /* At this point we've done that diagonalization and all that remains is
     * to certify the correctness and compute error bounds. */
    {
        arb_mat_t e;
        arb_struct *error_norms;
        arf_t x1, x2;
        arb_t z1, z2;
        arb_mat_t B;

        error_norms = _arb_vec_init(dim);

        arb_mat_init(B, dim, dim);
        arb_mat_init(e, dim, 1);

        arf_init(x1);
        arf_init(x2);
        arb_init(z1);
        arb_init(z2);

        for (j = 0; j < dim; j++)
        {
            arb_mat_set(B, A);
            for (k = 0; k < dim; k++)
            {
                arb_sub(B(k, k), B(k, k), D(j), prec);
            }
            for (k = 0; k < dim; k++)
            {
                arb_set(arb_mat_entry(e, k, 0), P(k, j));
            }
            arb_mat_frobenius_norm(z2, e, prec);
            arb_mat_mul(e, B, e, prec);
            arb_mat_frobenius_norm(error_norms+j, e, prec);

            /* And now z2 is an upper bound for the
             * error in the eigenvalue. */
            arb_div(z2, error_norms+j, z2, prec);
            arb_add_error(D(j), z2);
        }

        for (j = 0; j < dim; j++)
        {
            if (j == 0)
            {
                arb_sub(z1, D(j), D(1), prec);
            }
            else
            {
                arb_sub(z1, D(j), D(0), prec);
            }
            arb_get_abs_lbound_arf(x1, z1, prec);
            for (k = 1; k < dim; k++)
            {
                if (k == j)
                    continue;
                arb_sub(z1, D(j), D(k), prec);
                arb_get_abs_lbound_arf(x2, z1, prec);
                if (arf_cmp(x2, x1) < 0)
                {
                    arf_set(x1, x2);
                }
            }
            if (arf_is_zero(x1))
            {
                unique_eigenvalues = 0;
            }
            arb_div_arf(z1, error_norms+j, x1, prec);
            for (k = 0; k < dim; k++)
            {
                arb_add_error(P(k, j), z1);
            }
        }

        arb_mat_clear(B);
        arb_mat_clear(e);
        _arb_vec_clear(error_norms, dim);
        arf_clear(x1);
        arf_clear(x2);
        arb_clear(z1);
        arb_clear(z2);
    }

    _sort_decomposition(D, P);
    _standardize_column_norms(P, prec);
    _standardize_column_signs(P);

    if (unique_eigenvalues)
        return 0;
    else
        return 1;

#undef B
#undef D
#undef P
}
