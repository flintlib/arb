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

#include "fmpcb_poly.h"

#define CUTOFF 5

void
_fmpcb_poly_revert_series_newton(fmpcb_ptr Qinv, fmpcb_srcptr Q, long Qlen, long n, long prec)
{
    long i, k, a[FLINT_BITS];
    fmpcb_ptr T, U, V;

    if (n <= 2)
    {
        if (n >= 1)
            fmpcb_zero(Qinv);
        if (n == 2)
            fmpcb_inv(Qinv + 1, Q + 1, prec);
        return;
    }

    T = _fmpcb_vec_init(n);
    U = _fmpcb_vec_init(n);
    V = _fmpcb_vec_init(n);

    k = n;
    for (i = 1; (1L << i) < k; i++);
    a[i = 0] = k;
    while (k >= CUTOFF)
        a[++i] = (k = (k + 1) / 2);

    _fmpcb_poly_revert_series_lagrange(Qinv, Q, Qlen, k, prec);
    _fmpcb_vec_zero(Qinv + k, n - k);

    for (i--; i >= 0; i--)
    {
        k = a[i];
        _fmpcb_poly_compose_series(T, Q, FLINT_MIN(Qlen, k), Qinv, k, k, prec);
        _fmpcb_poly_derivative(U, T, k, prec); fmpcb_zero(U + k - 1);
        fmpcb_zero(T + 1);
        _fmpcb_poly_div_series(V, T, k, U, k, k, prec);
        _fmpcb_poly_derivative(T, Qinv, k, prec);
        _fmpcb_poly_mullow(U, V, k, T, k, k, prec);
        _fmpcb_vec_sub(Qinv, Qinv, U, k, prec);
    }

    _fmpcb_vec_clear(T, n);
    _fmpcb_vec_clear(U, n);
    _fmpcb_vec_clear(V, n);
}

void
fmpcb_poly_revert_series_newton(fmpcb_poly_t Qinv,
                                    const fmpcb_poly_t Q, long n, long prec)
{
    long Qlen = Q->length;

    if (Qlen < 2 || !fmpcb_is_zero(Q->coeffs)
                 || fmpcb_contains_zero(Q->coeffs + 1))
    {
        printf("Exception (fmpcb_poly_revert_series_newton). Input must \n"
               "have zero constant term and nonzero coefficient of x^1.\n");
        abort();
    }

    if (Qinv != Q)
    {
        fmpcb_poly_fit_length(Qinv, n);
        _fmpcb_poly_revert_series_newton(Qinv->coeffs, Q->coeffs, Qlen, n, prec);
    }
    else
    {
        fmpcb_poly_t t;
        fmpcb_poly_init2(t, n);
        _fmpcb_poly_revert_series_newton(t->coeffs, Q->coeffs, Qlen, n, prec);
        fmpcb_poly_swap(Qinv, t);
        fmpcb_poly_clear(t);
    }

    _fmpcb_poly_set_length(Qinv, n);
    _fmpcb_poly_normalise(Qinv);
}

