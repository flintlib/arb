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
#include "gamma.h"
#include "zeta.h"

static __inline__ void
_log_rfac_series(fmpcb_ptr t, const fmpcb_t x, long r, long len, long prec)
{
    fmpcb_struct f[2];
    long i, rflen;

    fmpcb_init(f);
    fmpcb_init(f + 1);
    fmpcb_set(f, x);
    fmpcb_one(f + 1);

    rflen = FLINT_MIN(len, r + 1);
    _fmpcb_poly_rfac_series_ui(t, f, FLINT_MIN(2, len), r, rflen, prec);
    _fmpcb_poly_log_series(t, t, rflen, len, prec);

    /* now get the right branch cut for the constant term
       TODO: make this a proper function */
    fmpcb_zero(t);
    for (i = 0; i < r; i++)
    {
        fmpcb_add_ui(f, x, i, prec);
        fmpcb_log(f, f, prec);
        fmpcb_add(t, t, f, prec);
    }

    fmpcb_clear(f);
    fmpcb_clear(f + 1);
}

void
_fmpcb_poly_lgamma_series(fmpcb_ptr res, fmpcb_srcptr h, long hlen, long len, long prec)
{
    int reflect;
    long r, n, wp;
    fmpcb_t zr;
    fmpcb_ptr t, u;

    hlen = FLINT_MIN(hlen, len);
    wp = prec + FLINT_BIT_COUNT(prec);

    t = _fmpcb_vec_init(len);
    u = _fmpcb_vec_init(len);
    fmpcb_init(zr);

    /* TODO: use real code at real numbers */
    if (0)
    {
    }
    else
    {
        /* otherwise use Stirling series */
        gamma_stirling_choose_param_fmpcb(&reflect, &r, &n, h, 0, 0, wp);
        fmpcb_add_ui(zr, h, r, wp);
        gamma_stirling_eval_fmpcb_series(u, zr, n, len, wp);

        if (r != 0)
        {
            _log_rfac_series(t, h, r, len, wp);
            _fmpcb_vec_sub(u, u, t, len, wp);
        }
    }

    /* compose with nonconstant part */
    fmpcb_zero(t);
    _fmpcb_vec_set(t + 1, h + 1, hlen - 1);
    _fmpcb_poly_compose_series(res, u, len, t, hlen, len, prec);

    fmpcb_clear(zr);
    _fmpcb_vec_clear(t, len);
    _fmpcb_vec_clear(u, len);
}

void
fmpcb_poly_lgamma_series(fmpcb_poly_t res, const fmpcb_poly_t f, long n, long prec)
{
    fmpcb_poly_fit_length(res, n);

    if (f->length == 0 || n == 0)
        _fmpcb_vec_indeterminate(res->coeffs, n);
    else
        _fmpcb_poly_lgamma_series(res->coeffs, f->coeffs, f->length, n, prec);

    _fmpcb_poly_set_length(res, n);
    _fmpcb_poly_normalise(res);
}

