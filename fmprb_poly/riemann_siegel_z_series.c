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
_fmprb_poly_riemann_siegel_z_series(fmprb_ptr res, fmprb_srcptr h, long hlen, long len, long prec)
{
    long i, alloc;
    fmprb_ptr t, u, v, w, q;

    hlen = FLINT_MIN(hlen, len);

    alloc = 5 * len;
    t = _fmprb_vec_init(alloc);
    u = t + len;
    v = u + len;
    w = v + len;
    q = w + len;

    /* (v + wi) = zeta(1/2 + i (s+x)) */
    {
        fmpcb_ptr sx, z;
        fmpcb_t a;
        long slen = FLINT_MIN(len, 2);

        z = _fmpcb_vec_init(len);
        sx = _fmpcb_vec_init(slen);
        fmpcb_init(a);

        fmpcb_one(a);
        fmpcb_one(sx);
        fmpcb_mul_2exp_si(sx, sx, -1);
        fmprb_set(fmpcb_imagref(sx), h);
        if (slen > 1)
            fmprb_one(fmpcb_imagref(sx + 1));

        _fmpcb_poly_zeta_series(z, sx, slen, a, 0, len, prec);

        for (i = 0; i < len; i++)
        {
            fmprb_set(v + i, fmpcb_realref(z + i));
            fmprb_set(w + i, fmpcb_imagref(z + i));
        }

        fmpcb_clear(a);
        _fmpcb_vec_clear(z, len);
        _fmpcb_vec_clear(sx, slen);
    }

    /* (t + ui) = exp(i theta(s+x)) */
    fmprb_set(u, h);
    if (len > 1)
        fmprb_one(u + 1);
    _fmprb_poly_riemann_siegel_theta_series(t, u, 2, len, prec);
    _fmprb_poly_sin_cos_series(u, t, t, len, len, prec);

    /* (t+ui)(v+wi) = (tv-uw) + (tw+uv)i */
    _fmprb_poly_mullow(q, t, len, v, len, len, prec);
    _fmprb_poly_mullow(t, u, len, w, len, len, prec);
    _fmprb_vec_sub(t, q, t, len, prec);

    /* compose with nonconstant part */
    fmprb_zero(u);
    _fmprb_vec_set(u + 1, h + 1, hlen - 1);
    _fmprb_poly_compose_series(res, t, len, u, hlen, len, prec);

    _fmprb_vec_clear(t, alloc);
}

void
fmprb_poly_riemann_siegel_z_series(fmprb_poly_t res, const fmprb_poly_t f, long n, long prec)
{
    if (n == 0)
    {
        fmprb_poly_zero(res);
        return;
    }

    fmprb_poly_fit_length(res, n);

    if (f->length == 0)
    {
        fmprb_t t;
        fmprb_init(t);
        _fmprb_poly_riemann_siegel_z_series(res->coeffs, t, 1, n, prec);
        fmprb_clear(t);
    }
    else
    {
        _fmprb_poly_riemann_siegel_z_series(res->coeffs, f->coeffs, f->length, n, prec);
    }

    _fmprb_poly_set_length(res, n);
    _fmprb_poly_normalise(res);
}

