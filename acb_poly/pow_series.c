/*=============================================================================

    This file is part of acb.

    acb is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    acb is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with acb; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2013 Fredrik Johansson

******************************************************************************/

#include "acb_poly.h"

void
_acb_poly_pow_series(acb_ptr h,
    acb_srcptr f, long flen,
    acb_srcptr g, long glen, long len, long prec)
{
    if (glen == 1)
    {
        _acb_poly_pow_acb_series(h, f, flen, g, len, prec);
        return;
    }

    /* f^g = exp(g * log(f)) */
    if (flen == 1)
    {
        acb_t t;
        acb_init(t);
        acb_log(t, f, prec);
        _acb_vec_scalar_mul(h, g, glen, t, prec);
        _acb_poly_exp_series(h, h, glen, len, prec);
        acb_clear(t);
    }
    else
    {
        acb_ptr t;
        t = _acb_vec_init(len);
        _acb_poly_log_series(t, f, flen, len, prec);
        _acb_poly_mullow(h, t, len, g, glen, len, prec);
        _acb_poly_exp_series(h, h, len, len, prec);
        _acb_vec_clear(t, len);

    }
}

void
acb_poly_pow_series(acb_poly_t h,
    const acb_poly_t f, const acb_poly_t g, long len, long prec)
{
    long flen, glen;

    flen = f->length;
    glen = g->length;

    flen = FLINT_MIN(flen, len);
    glen = FLINT_MIN(glen, len);

    if (len == 0)
    {
        acb_poly_zero(h);
        return;
    }

    if (glen == 0)
    {
        acb_poly_one(h);
        return;
    }

    if (flen == 0)
    {
        acb_poly_zero(h);
        return;
    }

    if (flen == 1 && glen == 1)
    {
        acb_poly_fit_length(h, 1);
        acb_pow(h->coeffs, f->coeffs, g->coeffs, prec);
        _acb_poly_set_length(h, 1);
        _acb_poly_normalise(h);
        return;
    }

    if (f == h || g == h)
    {
        acb_poly_t t;
        acb_poly_init2(t, len);
        _acb_poly_pow_series(t->coeffs,
            f->coeffs, flen, g->coeffs, glen, len, prec);
        _acb_poly_set_length(t, len);
        _acb_poly_normalise(t);
        acb_poly_swap(t, h);
        acb_poly_clear(t);
    }
    else
    {
        acb_poly_fit_length(h, len);
        _acb_poly_pow_series(h->coeffs,
            f->coeffs, flen, g->coeffs, glen, len, prec);
        _acb_poly_set_length(h, len);
        _acb_poly_normalise(h);
    }
}

