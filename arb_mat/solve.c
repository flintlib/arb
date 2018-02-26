/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"

int
arb_mat_solve_lu(arb_mat_t X, const arb_mat_t A, const arb_mat_t B, slong prec)
{
    int result;
    slong n, m, *perm;
    arb_mat_t LU;

    n = arb_mat_nrows(A);
    m = arb_mat_ncols(X);

    if (n == 0 || m == 0)
        return 1;

    perm = _perm_init(n);
    arb_mat_init(LU, n, n);

    result = arb_mat_lu(perm, LU, A, prec);

    if (result)
        arb_mat_solve_lu_precomp(X, perm, LU, B, prec);

    arb_mat_clear(LU);
    _perm_clear(perm);

    return result;
}

/*
 * Helper function to compute a lower bound of 1 - inf_norm(I - A*B).
 * Returns zero when this lower bound is zero.
 */
int _mag_err_complement(mag_t m,
    const arb_mat_t A, const arb_mat_t B, slong prec)
{
    slong n;
    arb_mat_t AB, E;
    mag_t err;

    n = arb_mat_nrows(A);

    mag_init(err);
    arb_mat_init(AB, n, n);
    arb_mat_init(E, n, n);

    arb_mat_mul(AB, A, B, prec);
    arb_mat_one(E);
    arb_mat_sub(E, E, AB, prec);
    arb_mat_bound_inf_norm(err, E);
    mag_one(m);
    mag_sub_lower(m, m, err);

    mag_clear(err);
    arb_mat_clear(AB);
    arb_mat_clear(E);

    return !mag_is_zero(m);
}

int arb_mat_solve_precond_precomp(arb_mat_t X, const arb_mat_t A,
    const arb_mat_t B, const arb_mat_t R, const arb_mat_t T, slong prec)
{
    int result;
    slong m, n;
    mag_t d;

    result = 0;

    n = arb_mat_nrows(A);
    m = arb_mat_ncols(X);

    if (n == 0 || m == 0)
        return 1;

    /* Use Theorem 10.2 of Rump in Acta Numerica 2010 */
    mag_init(d);
    if (_mag_err_complement(d, R, A, prec))
    {
        arb_mat_t C;

        arb_mat_init(C, n, m);
        {
            arb_mat_t B_prime, B_error;

            arb_mat_init(B_prime, n, m);
            arb_mat_init(B_error, n, m);

            arb_mat_mul(B_prime, A, T, prec);
            arb_mat_sub(B_error, B, B_prime, prec);
            arb_mat_mul(C, R, B_error, prec);

            arb_mat_clear(B_prime);
            arb_mat_clear(B_error);
        }

        /* Each column gets its own error bound. */
        arb_mat_set(X, T);
        {
            int i, j;
            mag_t e, err;

            mag_init(e);
            mag_init(err);

            for (j = 0; j < m; j++)
            {
                mag_zero(err);
                for (i = 0; i < n; i++)
                {
                    arb_get_mag(e, arb_mat_entry(C, i, j));
                    mag_max(err, err, e);
                }
                mag_div(err, err, d);
                for (i = 0; i < n; i++)
                {
                    arb_add_error_mag(arb_mat_entry(X, i, j), err);
                }
            }

            mag_clear(e);
            mag_clear(err);
        }

        arb_mat_clear(C);

        result = 1;
    }

    mag_clear(d);

    return result;
}

int
arb_mat_solve_precond(arb_mat_t X,
    const arb_mat_t A, const arb_mat_t B, slong prec)
{
    int result;
    slong m, n;
    arb_mat_t R, T;

    n = arb_mat_nrows(A);
    m = arb_mat_ncols(X);

    if (n == 0 || m == 0)
        return 1;

    arb_mat_init(R, n, n);
    arb_mat_init(T, n, m);
    {
        slong *perm;
        arb_mat_t I, LU;

        perm = _perm_init(n);
        arb_mat_init(I, n, n);
        arb_mat_init(LU, n, n);

        arb_mat_one(I);
        result = arb_mat_approx_lu(perm, LU, A, prec);
        if (result)
        {
            arb_mat_approx_solve_lu_precomp(R, perm, LU, I, prec);
            arb_mat_approx_solve_lu_precomp(T, perm, LU, B, prec);
        }

        _perm_clear(perm);
        arb_mat_clear(I);
        arb_mat_clear(LU);
    }

    if (result)
        result = arb_mat_solve_precond_precomp(X, A, B, R, T, prec);

    arb_mat_clear(R);
    arb_mat_clear(T);

    return result;
}

int
arb_mat_solve_c(arb_mat_t X, const arb_mat_t A, const arb_mat_t B, slong prec)
{
    int result;
    slong m, n;
    arb_mat_t I, R;

    n = arb_mat_nrows(A);
    m = arb_mat_ncols(X);

    if (n == 0 || m == 0)
        return 1;

    arb_mat_init(I, n, n);
    arb_mat_init(R, n, n);

    arb_mat_one(I);
    result = arb_mat_approx_solve(R, A, I, prec);
    if (result)
    {
        arb_mat_t RA, RB;

        arb_mat_init(RA, n, n);
        arb_mat_init(RB, n, m);

        arb_mat_mul(RA, R, A, prec);
        arb_mat_mul(RB, R, B, prec);

        result = arb_mat_solve_lu(X, RA, RB, prec);

        arb_mat_clear(RA);
        arb_mat_clear(RB);
    }

    arb_mat_clear(I);
    arb_mat_clear(R);

    return result;
}

int
arb_mat_solve_d(arb_mat_t X, const arb_mat_t A, const arb_mat_t B, slong prec)
{
    int result;
    slong m, n;
    arb_mat_t I, R;

    n = arb_mat_nrows(A);
    m = arb_mat_ncols(X);

    if (n == 0 || m == 0)
        return 1;

    arb_mat_init(I, n, n);
    arb_mat_init(R, n, n);

    arb_mat_one(I);
    result = arb_mat_approx_solve(R, A, I, prec);
    if (result)
    {
        arb_mat_t RA, RB, E;
        mag_t d;

        arb_mat_init(RA, n, n);
        arb_mat_init(RB, n, m);
        arb_mat_init(E, n, n);
        mag_init(d);

        arb_mat_mul(RA, R, A, prec);
        arb_mat_mul(RB, R, B, prec);
        arb_mat_sub(E, I, RA, prec);
        arb_mat_bound_inf_norm(d, E);

        if (mag_cmp_2exp_si(d, 0) < 0)
        {
            int i, j;
            mag_t e, err;

            mag_init(e);
            mag_init(err);

            mag_geom_series(d, d, 1);
            arb_mat_set(X, RB);

            for (j = 0; j < m; j++)
            {
                mag_zero(err);
                for (i = 0; i < n; i++)
                {
                    arb_get_mag(e, arb_mat_entry(RB, i, j));
                    mag_max(err, err, e);
                }
                mag_mul(err, err, d);
                for (i = 0; i < n; i++)
                {
                    arb_add_error_mag(arb_mat_entry(X, i, j), err);
                }
            }

            mag_clear(e);
            mag_clear(err);
        }
        else
        {
            result = arb_mat_solve_lu(X, RA, RB, prec);
        }

        arb_mat_clear(RA);
        arb_mat_clear(RB);
        arb_mat_clear(E);
        mag_clear(d);
    }

    arb_mat_clear(I);
    arb_mat_clear(R);

    return result;
}

int
arb_mat_solve(arb_mat_t X, const arb_mat_t A, const arb_mat_t B, slong prec)
{
    return arb_mat_solve_lu(X, A, B, prec);
}
