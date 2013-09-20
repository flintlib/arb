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
#include "fmpcb_poly.h"

void
_fmprb_poly_riemann_siegel_theta_series(fmprb_ptr res,
    fmprb_srcptr h, long hlen, long len, long prec)
{
    fmpcb_ptr s;
    fmprb_t u;
    long i;

    hlen = FLINT_MIN(hlen, len);

    s = _fmpcb_vec_init(len);

    fmprb_init(u);

    /* s = 1/4 + (1/2) i h */
    for (i = 0; i < hlen; i++)
        fmprb_mul_2exp_si(fmpcb_imagref(s + i), h + i, -1);

    fmprb_one(u);
    fmprb_mul_2exp_si(u, u, -2);
    fmprb_add(fmpcb_realref(s), fmpcb_realref(s), u, prec);

    /* log gamma */
    _fmpcb_poly_lgamma_series(s, s, hlen, len, prec);

    /* imaginary part */
    for (i = 0; i < len; i++)
        fmprb_set(res + i, fmpcb_imagref(s + i));

    /* subtract log(pi)/2 * h */
    fmprb_const_pi(u, prec);
    fmprb_log(u, u, prec);
    fmprb_mul_2exp_si(u, u, -1);
    fmprb_neg(u, u);
    _fmprb_vec_scalar_addmul(res, h, hlen, u, prec);

    _fmpcb_vec_clear(s, len);
    fmprb_clear(u);
}

void
fmprb_poly_riemann_siegel_theta_series(fmprb_poly_t res,
    const fmprb_poly_t f, long n, long prec)
{
    if (n == 0 || f->length == 0)
    {
        fmprb_poly_zero(res);
        return;
    }

    if (res == f)
    {
        fmprb_poly_t tmp;
        fmprb_poly_init2(tmp, n);
        _fmprb_poly_riemann_siegel_theta_series(tmp->coeffs,
            f->coeffs, f->length, n, prec);
        fmprb_poly_swap(res, tmp);
        fmprb_poly_clear(tmp);
    }
    else
    {
        fmprb_poly_fit_length(res, n);
        _fmprb_poly_riemann_siegel_theta_series(res->coeffs,
            f->coeffs, f->length, n, prec);
    }

    _fmprb_poly_set_length(res, n);
    _fmprb_poly_normalise(res);
}

