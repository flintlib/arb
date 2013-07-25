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

/* pointer to (x/Q)^i */
#define Ri(ii) (R + (n-1)*((ii)-1))

void
_fmprb_poly_revert_series_lagrange_fast(fmprb_ptr Qinv, fmprb_srcptr Q, long n, long prec)
{
    long i, j, k, m;
    fmprb_ptr R, S, T, tmp;
    fmprb_t t;

    if (n <= 2)
    {
        if (n >= 1)
            fmprb_zero(Qinv);
        if (n == 2)
            fmprb_ui_div(Qinv + 1, 1, Q + 1, prec);
        return;
    }

    m = n_sqrt(n);

    fmprb_init(t);
    R = _fmprb_vec_init((n - 1) * m);
    S = _fmprb_vec_init(n - 1);
    T = _fmprb_vec_init(n - 1);

    fmprb_zero(Qinv);
    fmprb_ui_div(Qinv + 1, 1, Q + 1, prec);

    _fmprb_poly_inv_series(Ri(1), Q + 1, n - 1, n - 1, prec);
    for (i = 2; i <= m; i++)
        _fmprb_poly_mullow(Ri(i), Ri((i + 1) / 2), n - 1, Ri(i / 2), n - 1, n - 1, prec);

    for (i = 2; i < m; i++)
        fmprb_div_ui(Qinv + i, Ri(i) + i - 1, i, prec);

    _fmprb_vec_set(S, Ri(m), n - 1);

    for (i = m; i < n; i += m)
    {
        fmprb_div_ui(Qinv + i, S + i - 1, i, prec);

        for (j = 1; j < m && i + j < n; j++)
        {
            fmprb_mul(t, S + 0, Ri(j) + i + j - 1, prec);
            for (k = 1; k <= i + j - 1; k++)
                fmprb_addmul(t, S + k, Ri(j) + i + j - 1 - k, prec);
            fmprb_div_ui(Qinv + i + j, t, i + j, prec);
        }

        if (i + 1 < n)
        {
            _fmprb_poly_mullow(T, S, n - 1, Ri(m), n - 1, n - 1, prec);
            tmp = S; S = T; T = tmp;
        }
    }

    fmprb_clear(t);
    _fmprb_vec_clear(R, (n - 1) * m);
    _fmprb_vec_clear(S, n - 1);
    _fmprb_vec_clear(T, n - 1);
}

void
fmprb_poly_revert_series_lagrange_fast(fmprb_poly_t Qinv,
                                    const fmprb_poly_t Q, long n, long prec)
{
    fmprb_ptr Qcopy;
    int Qalloc;
    long Qlen = Q->length;

    if (Q->length < 2 || !fmprb_is_zero(Q->coeffs)
                      || fmprb_contains_zero(Q->coeffs + 1))
    {
        printf("Exception (fmprb_poly_revert_series_lagrange_fast). Input \n"
               "must have zero constant term and nonzero coefficient of x^1.\n");
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
        _fmprb_poly_revert_series_lagrange_fast(Qinv->coeffs, Qcopy, n, prec);
    }
    else
    {
        fmprb_poly_t t;
        fmprb_poly_init2(t, n);
        _fmprb_poly_revert_series_lagrange_fast(t->coeffs, Qcopy, n, prec);
        fmprb_poly_swap(Qinv, t);
        fmprb_poly_clear(t);
    }

    _fmprb_poly_set_length(Qinv, n);
    _fmprb_poly_normalise(Qinv);

    if (Qalloc)
        flint_free(Qcopy);
}

