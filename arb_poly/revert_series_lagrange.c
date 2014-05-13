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

#include "arb_poly.h"

void
_arb_poly_revert_series_lagrange(arb_ptr Qinv,
    arb_srcptr Q, long Qlen, long n, long prec)
{
    long i;
    arb_ptr R, S, T, tmp;

    if (n <= 2)
    {
        if (n >= 1)
            arb_zero(Qinv);
        if (n == 2)
            arb_inv(Qinv + 1, Q + 1, prec);
        return;
    }

    R = _arb_vec_init(n - 1);
    S = _arb_vec_init(n - 1);
    T = _arb_vec_init(n - 1);

    arb_zero(Qinv);
    arb_inv(Qinv + 1, Q + 1, prec);

    _arb_poly_inv_series(R, Q + 1, FLINT_MIN(Qlen, n) - 1, n - 1, prec);
    _arb_vec_set(S, R, n - 1);

    for (i = 2; i < n; i++)
    {
        _arb_poly_mullow(T, S, n - 1, R, n - 1, n - 1, prec);
        arb_div_ui(Qinv + i, T + i - 1, i, prec);
        tmp = S; S = T; T = tmp;
    }

    _arb_vec_clear(R, n - 1);
    _arb_vec_clear(S, n - 1);
    _arb_vec_clear(T, n - 1);
}

void
arb_poly_revert_series_lagrange(arb_poly_t Qinv,
                                    const arb_poly_t Q, long n, long prec)
{
    long Qlen = Q->length;

    if (Qlen < 2 || !arb_is_zero(Q->coeffs)
                 || arb_contains_zero(Q->coeffs + 1))
    {
        printf("Exception (arb_poly_revert_series_lagrange). Input must \n"
               "have zero constant term and nonzero coefficient of x^1.\n");
        abort();
    }

    if (Qinv != Q)
    {
        arb_poly_fit_length(Qinv, n);
        _arb_poly_revert_series_lagrange(Qinv->coeffs, Q->coeffs, Qlen, n, prec);
    }
    else
    {
        arb_poly_t t;
        arb_poly_init2(t, n);
        _arb_poly_revert_series_lagrange(t->coeffs, Q->coeffs, Qlen, n, prec);
        arb_poly_swap(Qinv, t);
        arb_poly_clear(t);
    }

    _arb_poly_set_length(Qinv, n);
    _arb_poly_normalise(Qinv);
}

