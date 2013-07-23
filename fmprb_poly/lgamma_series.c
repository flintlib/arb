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
#include "gamma.h"
#include "zeta.h"

static __inline__ int
fmprb_is_int(const fmprb_t x)
{
    return fmprb_is_zero(x) ||
        (fmprb_is_exact(x) &&
                 fmpz_sgn(fmpr_expref(fmprb_midref(x))) >= 0);
}

static __inline__ void
_fmprb_vec_indeterminate(fmprb_ptr vec, long len)
{
    long i;
    for (i = 0; i < len; i++)
    {
        fmpr_nan(fmprb_midref(vec + i));
        fmpr_pos_inf(fmprb_radref(vec + i));
    }
}

static __inline__ void
_log_rfac_series(fmprb_ptr t, const fmprb_t x, long r, long len, long prec)
{
    fmprb_struct f[2];
    long rflen;

    fmprb_init(f);
    fmprb_init(f + 1);
    fmprb_set(f, x);
    fmprb_one(f + 1);

    rflen = FLINT_MIN(len, r + 1);
    _fmprb_poly_rfac_series_ui(t, f, FLINT_MIN(2, len), r, rflen, prec);
    _fmprb_poly_log_series(t, t, rflen, len, prec);

    fmprb_clear(f);
    fmprb_clear(f + 1);
}

void
_fmprb_poly_lgamma_series(fmprb_ptr res, fmprb_ptr h, long hlen, long len, long prec)
{
    int reflect;
    long i, r, n, wp;
    fmprb_t zr;
    fmprb_ptr t, u;

    hlen = FLINT_MIN(hlen, len);
    wp = prec + FLINT_BIT_COUNT(prec);

    t = _fmprb_vec_init(len);
    u = _fmprb_vec_init(len);
    fmprb_init(zr);

    /* use zeta values at small integers */
    if (fmprb_is_int(h) && (fmpr_cmpabs_ui(fmprb_midref(h), prec / 2) < 0))
    {
        r = fmpr_get_si(fmprb_midref(h), FMPR_RND_DOWN);

        if (r <= 0)
        {
            _fmprb_vec_indeterminate(res, len);
        }
        else
        {
            fmprb_zero(u);
            if (len > 1) fmprb_const_euler(u + 1, wp);
            if (len > 2) zeta_ui_vec(u + 2, 2, len - 2, wp);
            for (i = 2; i < len; i++)
                fmprb_div_ui(u + i, u + i, i, wp);
            for (i = 1; i < len; i += 2)
                fmprb_neg(u + i, u + i);

            if (r != 1)
            {
                fmprb_one(zr);
                _log_rfac_series(t, zr, r - 1, len, wp);
                _fmprb_vec_add(u, u, t, len, wp);
            }
        }
    }
    else
    {
        /* otherwise use Stirling series */
        gamma_stirling_choose_param_fmprb(&reflect, &r, &n, h, 0, 0, wp);
        fmprb_add_ui(zr, h, r, wp);
        gamma_stirling_eval_fmprb_series(u, zr, n, len, wp);

        if (r != 0)
        {
            _log_rfac_series(t, h, r, len, wp);
            _fmprb_vec_sub(u, u, t, len, wp);
        }
    }

    /* compose with nonconstant part */
    fmprb_zero(t);
    _fmprb_vec_set(t + 1, h + 1, hlen - 1);
    _fmprb_poly_compose_series(res, u, len, t, hlen, len, prec);

    fmprb_clear(zr);
    _fmprb_vec_clear(t, len);
    _fmprb_vec_clear(u, len);
}

void
fmprb_poly_lgamma_series(fmprb_poly_t res, const fmprb_poly_t f, long n, long prec)
{
    if (f->length == 0 || n == 0)
    {
        printf("fmprb_poly_lgamma_series: require n > 0 and nonzero input\n");
        abort();
    }

    fmprb_poly_fit_length(res, n);
    _fmprb_poly_lgamma_series(res->coeffs, f->coeffs, f->length, n, prec);
    _fmprb_poly_set_length(res, n);
    _fmprb_poly_normalise(res);
}

