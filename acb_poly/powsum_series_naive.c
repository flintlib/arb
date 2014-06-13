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

    Copyright (C) 2014 Fredrik Johansson

******************************************************************************/

#include "acb_poly.h"

void _acb_poly_acb_invpow_cpx(acb_ptr res, const acb_t N, const acb_t c, long trunc, long prec)
{
    long i;
    acb_t logN;

    acb_init(logN);
    acb_log(logN, N, prec);
    acb_mul(res + 0, logN, c, prec);
    acb_neg(res + 0, res + 0);
    acb_exp(res + 0, res + 0, prec);

    for (i = 1; i < trunc; i++)
    {
        acb_mul(res + i, res + i - 1, logN, prec);
        acb_div_si(res + i, res + i, -i, prec);
    }

    acb_clear(logN);
}

void
_acb_poly_powsum_series_naive(acb_ptr z,
    const acb_t s, const acb_t a, long n, long len, long prec)
{
    long k;
    acb_ptr t;
    acb_t c;

    t = _acb_vec_init(len);
    acb_init(c);

    _acb_vec_zero(z, len);

    for (k = 0; k < n; k++)
    {
        acb_add_ui(c, a, k, prec);
        _acb_poly_acb_invpow_cpx(t, c, s, len, prec);
        _acb_vec_add(z, z, t, len, prec);
    }

    _acb_vec_clear(t, len);
    acb_clear(c);
}

