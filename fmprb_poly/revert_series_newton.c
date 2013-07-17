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

    Copyright (C) 2013 Fredrik Johansson

******************************************************************************/

#include "fmprb_poly.h"

#define CUTOFF 5

void
_fmprb_poly_revert_series_newton(fmprb_ptr Qinv, fmprb_srcptr Q, long n, long prec)
{
    long i, k, a[FLINT_BITS];
    fmprb_ptr T, U, V;

    if (n <= 2)
    {
        if (n >= 1)
            fmprb_zero(Qinv);
        if (n == 2)
            fmprb_ui_div(Qinv + 1, 1, Q + 1, prec);
        return;
    }

    T = _fmprb_vec_init(n);
    U = _fmprb_vec_init(n);
    V = _fmprb_vec_init(n);

    k = n;
    for (i = 1; (1L << i) < k; i++);
    a[i = 0] = k;
    while (k >= CUTOFF)
        a[++i] = (k = (k + 1) / 2);

    _fmprb_poly_revert_series_lagrange(Qinv, Q, k, prec);
    _fmprb_vec_zero(Qinv + k, n - k);

    for (i--; i >= 0; i--)
    {
        k = a[i];
        _fmprb_poly_compose_series(T, Q, k, Qinv, k, k, prec);
        _fmprb_poly_derivative(U, T, k, prec); fmprb_zero(U + k - 1);
        fmprb_zero(T + 1);
        _fmprb_poly_div_series(V, T, k, U, k, k, prec);
        _fmprb_poly_derivative(T, Qinv, k, prec);
        _fmprb_poly_mullow(U, V, k, T, k, k, prec);
        _fmprb_vec_sub(Qinv, Qinv, U, k, prec);
    }

    _fmprb_vec_clear(T, n);
    _fmprb_vec_clear(U, n);
    _fmprb_vec_clear(V, n);
}

void
fmprb_poly_revert_series_newton(fmprb_poly_t Qinv,
                                    const fmprb_poly_t Q, long n, long prec)
{
    fmprb_ptr Qcopy;
    int Qalloc;
    long Qlen = Q->length;

    if (Q->length < 2 || !fmprb_is_zero(Q->coeffs)
                      || fmprb_contains_zero(Q->coeffs + 1))
    {
        printf("Exception (fmprb_poly_revert_series_newton). Input must \n"
               "have zero constant term and nonzero coefficient of x^1.\n");
        abort();
    }

    if (n < 2)
    {
        fmprb_poly_zero(Qinv);
        return;
    }

    if (Qlen >= n)
    {
        Qcopy = Q->coeffs;
        Qalloc = 0;
    }
    else
    {
        long i;
        Qcopy = _fmprb_vec_init(n);
        for (i = 0; i < Qlen; i++)
            Qcopy[i] = Q->coeffs[i];
        Qalloc = 1;
    }

    if (Qinv != Q)
    {
        fmprb_poly_fit_length(Qinv, n);
        _fmprb_poly_revert_series_newton(Qinv->coeffs, Qcopy, n, prec);
    }
    else
    {
        fmprb_poly_t t;
        fmprb_poly_init2(t, n);
        _fmprb_poly_revert_series_newton(t->coeffs, Qcopy, n, prec);
        fmprb_poly_swap(Qinv, t);
        fmprb_poly_clear(t);
    }

    _fmprb_poly_set_length(Qinv, n);
    _fmprb_poly_normalise(Qinv);

    if (Qalloc)
        flint_free(Qcopy);
}

