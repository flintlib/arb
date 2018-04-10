/*
    Copyright (C) 2018 arbguest

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

static void
_acb_approx_mul(acb_t res, const acb_t x, const acb_t y, slong prec)
{
    arf_complex_mul(arb_midref(acb_realref(res)), arb_midref(acb_imagref(res)),
        arb_midref(acb_realref(x)), arb_midref(acb_imagref(x)), 
        arb_midref(acb_realref(y)), arb_midref(acb_imagref(y)), prec, ARB_RND);
}

static void
_acb_approx_submul(acb_t z, const acb_t x, const acb_t y, acb_t t, slong prec)
{
    _acb_approx_mul(t, x, y, prec);

    arf_sub(arb_midref(acb_realref(z)),
        arb_midref(acb_realref(z)), 
        arb_midref(acb_realref(t)), prec, ARB_RND);
    arf_sub(arb_midref(acb_imagref(z)),
        arb_midref(acb_imagref(z)), 
        arb_midref(acb_imagref(t)), prec, ARB_RND);
}

/* note: the tmp variable t should have zero radius */
static void
_acb_approx_div(acb_t z, const acb_t x, const acb_t y, acb_t t, slong prec)
{
    arf_set(arb_midref(acb_realref(t)), arb_midref(acb_realref(y)));
    arf_set(arb_midref(acb_imagref(t)), arb_midref(acb_imagref(y)));

    acb_inv(t, t, prec);

    mag_zero(arb_radref(acb_realref(t)));
    mag_zero(arb_radref(acb_imagref(t)));

    _acb_approx_mul(z, x, t, prec);
}

void
acb_mat_approx_solve_lu_precomp(acb_mat_t X, const slong * perm,
    const acb_mat_t A, const acb_mat_t B, slong prec)
{
    acb_t t;
    slong i, j, c, n, m;

    n = acb_mat_nrows(X);
    m = acb_mat_ncols(X);

    if (X == B)
    {
        acb_ptr tmp = flint_malloc(sizeof(acb_struct) * n);

        for (c = 0; c < m; c++)
        {
            for (i = 0; i < n; i++)
                tmp[i] = B->rows[perm[i]][c];
            for (i = 0; i < n; i++)
                X->rows[i][c] = tmp[i];
        }

        flint_free(tmp);
    }
    else
    {
        for (c = 0; c < m; c++)
        {
            for (i = 0; i < n; i++)
            {
                acb_set(acb_mat_entry(X, i, c),
                    acb_mat_entry(B, perm[i], c));
            }
        }
    }

    acb_init(t);
    acb_mat_get_mid(X, X);

    for (c = 0; c < m; c++)
    {
        /* solve Ly = b */
        for (i = 1; i < n; i++)
        {
            for (j = 0; j < i; j++)
            {
                _acb_approx_submul(acb_mat_entry(X, i, c),
                    acb_mat_entry(A, i, j), acb_mat_entry(X, j, c), t, prec);
            }
        }

        /* solve Ux = y */
        for (i = n - 1; i >= 0; i--)
        {
            for (j = i + 1; j < n; j++)
            {
                _acb_approx_submul(acb_mat_entry(X, i, c),
                    acb_mat_entry(A, i, j), acb_mat_entry(X, j, c), t, prec);
            }

            _acb_approx_div(acb_mat_entry(X, i, c), acb_mat_entry(X, i, c),
                acb_mat_entry(A, i, i), t, prec);
        }
    }

    acb_clear(t);
}

