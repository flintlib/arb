/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

static void
_acb_mat_det_cofactor_2x2(acb_t t, const acb_mat_t A, slong prec)
{
    acb_mul   (t, acb_mat_entry(A, 0, 0), acb_mat_entry(A, 1, 1), prec);
    acb_submul(t, acb_mat_entry(A, 0, 1), acb_mat_entry(A, 1, 0), prec);
}

static void
_acb_mat_det_cofactor_3x3(acb_t t, const acb_mat_t A, slong prec)
{
    acb_t a;
    acb_init(a);

    acb_mul   (a, acb_mat_entry(A, 1, 0), acb_mat_entry(A, 2, 1), prec);
    acb_submul(a, acb_mat_entry(A, 1, 1), acb_mat_entry(A, 2, 0), prec);
    acb_mul   (t, a, acb_mat_entry(A, 0, 2), prec);

    acb_mul   (a, acb_mat_entry(A, 1, 2), acb_mat_entry(A, 2, 0), prec);
    acb_submul(a, acb_mat_entry(A, 1, 0), acb_mat_entry(A, 2, 2), prec);
    acb_addmul(t, a, acb_mat_entry(A, 0, 1), prec);

    acb_mul   (a, acb_mat_entry(A, 1, 1), acb_mat_entry(A, 2, 2), prec);
    acb_submul(a, acb_mat_entry(A, 1, 2), acb_mat_entry(A, 2, 1), prec);
    acb_addmul(t, a, acb_mat_entry(A, 0, 0), prec);

    acb_clear(a);
}

int
acb_mat_is_finite(const acb_mat_t A)
{
    slong i, j, n, m;

    n = acb_mat_nrows(A);
    m = acb_mat_ncols(A);

    for (i = 0; i < n; i++)
        for (j = 0; j < m; j++)
            if (!acb_is_finite(acb_mat_entry(A, i, j)))
                return 0;

    return 1;
}

int
acb_mat_is_triu(const acb_mat_t A)
{
    slong i, j, n, m;

    n = acb_mat_nrows(A);
    m = acb_mat_ncols(A);

    for (i = 1; i < n; i++)
        for (j = 0; j < FLINT_MIN(i, m); j++)
            if (!acb_is_zero(acb_mat_entry(A, i, j)))
                return 0;

    return 1;
}

int
acb_mat_is_tril(const acb_mat_t A)
{
    slong i, j, n, m;

    n = acb_mat_nrows(A);
    m = acb_mat_ncols(A);

    for (i = 0; i < n; i++)
        for (j = i + 1; j < m; j++)
            if (!acb_is_zero(acb_mat_entry(A, i, j)))
                return 0;

    return 1;
}

void
acb_mat_diag_prod(acb_t res, const acb_mat_t A, slong a, slong b, slong prec)
{
    if (b - a == 0)
        acb_one(res);
    else if (b - a == 1)
        acb_set_round(res, acb_mat_entry(A, a, a), prec);
    else if (b - a == 2)
        acb_mul(res, acb_mat_entry(A, a, a), acb_mat_entry(A, a + 1, a + 1), prec);
    else if (b - a == 3)
    {
        acb_mul(res, acb_mat_entry(A, a, a), acb_mat_entry(A, a + 1, a + 1), prec);
        acb_mul(res, res, acb_mat_entry(A, a + 2, a + 2), prec);
    }
    else
    {
        acb_t t;
        acb_init(t);
        acb_mat_diag_prod(t, A, a, a + (b - a) / 2, prec);
        acb_mat_diag_prod(res, A, a + (b - a) / 2, b, prec);
        acb_mul(res, res, t, prec);
        acb_clear(t);
    }
}

void
acb_mat_det(acb_t det, const acb_mat_t A, slong prec)
{
    slong n;

    if (!acb_mat_is_square(A))
    {
        flint_printf("acb_mat_det: a square matrix is required!\n");
        flint_abort();
    }

    n = acb_mat_nrows(A);

    if (n == 0)
    {
        acb_one(det);
    }
    else if (n == 1)
    {
        acb_set_round(det, acb_mat_entry(A, 0, 0), prec);
    }
    else if (n == 2)
    {
        _acb_mat_det_cofactor_2x2(det, A, prec);
    }
    else if (!acb_mat_is_finite(A))
    {
        acb_indeterminate(det);
    }
    else if (acb_mat_is_tril(A) || acb_mat_is_triu(A))
    {
        acb_mat_diag_prod(det, A, 0, n, prec);
    }
    else if (n == 3)
    {
        _acb_mat_det_cofactor_3x3(det, A, prec);
        /* note: 4x4 performs worse than LU */
    }
    else
    {
        if (n <= 14 || prec > 10.0 * n)
            acb_mat_det_lu(det, A, prec);
        else
            acb_mat_det_precond(det, A, prec);
    }
}

