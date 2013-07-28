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
    fmprb_t pi, u, v;
    fmpz_t pi_mult;
    long i, rflen, argprec;

    fmpcb_init(f);
    fmpcb_init(f + 1);
    fmprb_init(u);
    fmprb_init(pi);
    fmprb_init(v);
    fmpz_init(pi_mult);

    fmpcb_set(f, x);
    fmpcb_one(f + 1);

    rflen = FLINT_MIN(len, r + 1);
    _fmpcb_poly_rfac_series_ui(t, f, FLINT_MIN(2, len), r, rflen, prec);
    _fmpcb_poly_log_series(t, t, rflen, len, prec);

    /* now get the right branch cut for the constant term
       TODO: make this a proper function */
    argprec = FLINT_MIN(prec, 40);

    fmprb_zero(u);
    for (i = 0; i < r; i++)
    {
        fmpcb_add_ui(f, x, i, argprec);
        fmpcb_arg(v, f, argprec);
        fmprb_add(u, u, v, argprec);
    }

    if (argprec == prec)
    {
        fmprb_set(fmpcb_imagref(t), u);
    }
    else
    {
        fmprb_sub(v, u, fmpcb_imagref(t), argprec);
        fmprb_const_pi(pi, argprec);
        fmprb_div(v, v, pi, argprec);

        if (fmprb_get_unique_fmpz(pi_mult, v))
        {
            fmprb_const_pi(v, prec);
            fmprb_mul_fmpz(v, v, pi_mult, prec);
            fmprb_add(fmpcb_imagref(t), fmpcb_imagref(t), v, prec);
        }
        else
        {
            fmprb_zero(u);
            for (i = 0; i < r; i++)
            {
                fmpcb_add_ui(f, x, i, prec);
                fmpcb_arg(v, f, prec);
                fmprb_add(u, u, v, prec);
            }
            fmprb_set(fmpcb_imagref(t), u);
        }
    }

    fmpcb_clear(f);
    fmpcb_clear(f + 1);
    fmprb_clear(u);
    fmprb_clear(v);
    fmprb_clear(pi);
    fmpz_clear(pi_mult);
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

