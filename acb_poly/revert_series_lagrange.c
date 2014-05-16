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

void
_acb_poly_revert_series_lagrange(acb_ptr Qinv,
    acb_srcptr Q, long Qlen, long n, long prec)
{
    long i;
    acb_ptr R, S, T, tmp;

    if (n <= 2)
    {
        if (n >= 1)
            acb_zero(Qinv);
        if (n == 2)
            acb_inv(Qinv + 1, Q + 1, prec);
        return;
    }

    R = _acb_vec_init(n - 1);
    S = _acb_vec_init(n - 1);
    T = _acb_vec_init(n - 1);

    acb_zero(Qinv);
    acb_inv(Qinv + 1, Q + 1, prec);

    _acb_poly_inv_series(R, Q + 1, FLINT_MIN(Qlen, n) - 1, n - 1, prec);
    _acb_vec_set(S, R, n - 1);

    for (i = 2; i < n; i++)
    {
        _acb_poly_mullow(T, S, n - 1, R, n - 1, n - 1, prec);
        acb_div_ui(Qinv + i, T + i - 1, i, prec);
        tmp = S; S = T; T = tmp;
    }

    _acb_vec_clear(R, n - 1);
    _acb_vec_clear(S, n - 1);
    _acb_vec_clear(T, n - 1);
}

void
acb_poly_revert_series_lagrange(acb_poly_t Qinv,
                                    const acb_poly_t Q, long n, long prec)
{
    long Qlen = Q->length;

    if (Qlen < 2 || !acb_is_zero(Q->coeffs)
                 || acb_contains_zero(Q->coeffs + 1))
    {
        printf("Exception (acb_poly_revert_series_lagrange). Input must \n"
               "have zero constant term and nonzero coefficient of x^1.\n");
        abort();
    }

    if (Qinv != Q)
    {
        acb_poly_fit_length(Qinv, n);
        _acb_poly_revert_series_lagrange(Qinv->coeffs, Q->coeffs, Qlen, n, prec);
    }
    else
    {
        acb_poly_t t;
        acb_poly_init2(t, n);
        _acb_poly_revert_series_lagrange(t->coeffs, Q->coeffs, Qlen, n, prec);
        acb_poly_swap(Qinv, t);
        acb_poly_clear(t);
    }

    _acb_poly_set_length(Qinv, n);
    _acb_poly_normalise(Qinv);
}

