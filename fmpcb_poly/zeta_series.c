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

/* series of c^(d+x) */
static __inline__ void
_fmpcb_poly_pow_cpx(fmpcb_ptr res, const fmpcb_t c, const fmpcb_t d, long trunc, long prec)
{
    long i;
    fmpcb_t logc;

    fmpcb_init(logc);
    fmpcb_log(logc, c, prec);
    fmpcb_mul(res + 0, logc, d, prec);
    fmpcb_exp(res + 0, res + 0, prec);

    for (i = 1; i < trunc; i++)
    {
        fmpcb_mul(res + i, res + i - 1, logc, prec);
        fmpcb_div_ui(res + i, res + i, i, prec);
    }

    fmpcb_clear(logc);
}

void
_fmpcb_poly_zeta_series(fmpcb_ptr res, fmpcb_srcptr h, long hlen, const fmpcb_t a, int deflate, long len, long prec)
{
    long i;
    fmpcb_ptr t, u;

    hlen = FLINT_MIN(hlen, len);

    t = _fmpcb_vec_init(len);
    u = _fmpcb_vec_init(len);

    /* use reflection formula */
    if (fmpr_sgn(fmprb_midref(fmpcb_realref(h))) < 0 && fmpcb_is_one(a))
    {
        /* zeta(s) = (2*pi)**s * sin(pi*s/2) / pi * gamma(1-s) * zeta(1-s) */
        fmpcb_t pi;
        fmpcb_ptr f, s1, s2, s3, s4;

        fmpcb_init(pi);
        f = _fmpcb_vec_init(2);
        s1 = _fmpcb_vec_init(len);
        s2 = _fmpcb_vec_init(len);
        s3 = _fmpcb_vec_init(len);
        s4 = _fmpcb_vec_init(len);

        fmpcb_const_pi(pi, prec);

        /* s1 = (2*pi)**s */
        fmpcb_mul_2exp_si(pi, pi, 1);
        _fmpcb_poly_pow_cpx(s1, pi, h, len, prec);
        fmpcb_mul_2exp_si(pi, pi, -1);

        /* s2 = sin(pi*s/2) / pi */
        fmpcb_mul_2exp_si(pi, pi, -1);
        fmpcb_mul(f, pi, h, prec);
        fmpcb_set(f + 1, pi);
        fmpcb_mul_2exp_si(pi, pi, 1);
        _fmpcb_poly_sin_series(s2, f, 2, len, prec);
        _fmpcb_vec_scalar_div(s2, s2, len, pi, prec);

        /* s3 = gamma(1-s) */
        fmpcb_sub_ui(f, h, 1, prec);
        fmpcb_neg(f, f);
        fmpcb_set_si(f + 1, -1);
        _fmpcb_poly_gamma_series(s3, f, 2, len, prec);

        /* s4 = zeta(1-s) */
        fmpcb_sub_ui(f, h, 1, prec);
        fmpcb_neg(f, f);
        zeta_series(s4, f, a, 0, len, prec);
        for (i = 1; i < len; i += 2)
            fmpcb_neg(s4 + i, s4 + i);

        _fmpcb_poly_mullow(u, s1, len, s2, len, len, prec);
        _fmpcb_poly_mullow(s1, s3, len, s4, len, len, prec);
        _fmpcb_poly_mullow(t, u, len, s1, len, len, prec);

        /* add 1/(1-(s+t)) = 1/(1-s) + t/(1-s)^2 + ... */
        if (deflate)
        {
            fmpcb_sub_ui(u, h, 1, prec);
            fmpcb_neg(u, u);
            fmpcb_inv(u, u, prec);
            for (i = 1; i < len; i++)
                fmpcb_mul(u + i, u + i - 1, u, prec);
            _fmpcb_vec_add(t, t, u, len, prec);
        }

        fmpcb_clear(pi);
        _fmpcb_vec_clear(f, 2);
        _fmpcb_vec_clear(s1, len);
        _fmpcb_vec_clear(s2, len);
        _fmpcb_vec_clear(s3, len);
        _fmpcb_vec_clear(s4, len);
    }
    else
    {
        zeta_series(t, h, a, deflate, len, prec);
    }

    /* compose with nonconstant part */
    fmpcb_zero(u);
    _fmpcb_vec_set(u + 1, h + 1, hlen - 1);
    _fmpcb_poly_compose_series(res, t, len, u, hlen, len, prec);

    _fmpcb_vec_clear(t, len);
    _fmpcb_vec_clear(u, len);
}

void
fmpcb_poly_zeta_series(fmpcb_poly_t res, const fmpcb_poly_t f, const fmpcb_t a, int deflate, long n, long prec)
{
    if (n == 0)
    {
        fmpcb_poly_zero(res);
        return;
    }

    fmpcb_poly_fit_length(res, n);

    if (f->length == 0)
    {
        fmpcb_t t;
        fmpcb_init(t);
        _fmpcb_poly_zeta_series(res->coeffs, t, 1, a, deflate, n, prec);
        fmpcb_clear(t);
    }
    else
    {
        _fmpcb_poly_zeta_series(res->coeffs, f->coeffs, f->length, a, deflate, n, prec);
    }

    _fmpcb_poly_set_length(res, n);
    _fmpcb_poly_normalise(res);
}

