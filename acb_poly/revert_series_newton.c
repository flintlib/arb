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

#include "acb_poly.h"

#define CUTOFF 5

void
_acb_poly_revert_series_newton(acb_ptr Qinv, acb_srcptr Q, long Qlen, long n, long prec)
{
    long i, k, a[FLINT_BITS];
    acb_ptr T, U, V;

    if (n <= 2)
    {
        if (n >= 1)
            acb_zero(Qinv);
        if (n == 2)
            acb_inv(Qinv + 1, Q + 1, prec);
        return;
    }

    T = _acb_vec_init(n);
    U = _acb_vec_init(n);
    V = _acb_vec_init(n);

    k = n;
    for (i = 1; (1L << i) < k; i++);
    a[i = 0] = k;
    while (k >= CUTOFF)
        a[++i] = (k = (k + 1) / 2);

    _acb_poly_revert_series_lagrange(Qinv, Q, Qlen, k, prec);
    _acb_vec_zero(Qinv + k, n - k);

    for (i--; i >= 0; i--)
    {
        k = a[i];
        _acb_poly_compose_series(T, Q, FLINT_MIN(Qlen, k), Qinv, k, k, prec);
        _acb_poly_derivative(U, T, k, prec); acb_zero(U + k - 1);
        acb_zero(T + 1);
        _acb_poly_div_series(V, T, k, U, k, k, prec);
        _acb_poly_derivative(T, Qinv, k, prec);
        _acb_poly_mullow(U, V, k, T, k, k, prec);
        _acb_vec_sub(Qinv, Qinv, U, k, prec);
    }

    _acb_vec_clear(T, n);
    _acb_vec_clear(U, n);
    _acb_vec_clear(V, n);
}

void
acb_poly_revert_series_newton(acb_poly_t Qinv,
                                    const acb_poly_t Q, long n, long prec)
{
    long Qlen = Q->length;

    if (Qlen < 2 || !acb_is_zero(Q->coeffs)
                 || acb_contains_zero(Q->coeffs + 1))
    {
        printf("Exception (acb_poly_revert_series_newton). Input must \n"
               "have zero constant term and nonzero coefficient of x^1.\n");
        abort();
    }

    if (Qinv != Q)
    {
        acb_poly_fit_length(Qinv, n);
        _acb_poly_revert_series_newton(Qinv->coeffs, Q->coeffs, Qlen, n, prec);
    }
    else
    {
        acb_poly_t t;
        acb_poly_init2(t, n);
        _acb_poly_revert_series_newton(t->coeffs, Q->coeffs, Qlen, n, prec);
        acb_poly_swap(Qinv, t);
        acb_poly_clear(t);
    }

    _acb_poly_set_length(Qinv, n);
    _acb_poly_normalise(Qinv);
}

