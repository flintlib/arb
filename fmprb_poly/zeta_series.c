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

static __inline__ void
_fmprb_vec_printd(fmprb_srcptr vec, long len, long digits)
{
    long i;
    for (i = 0; i < len; i++)
        fmprb_printd(vec + i, digits), printf("\n");
}


/* series of c^(d+x) */
static __inline__ void
_fmprb_poly_pow_cpx(fmprb_ptr res, const fmprb_t c, const fmprb_t d, long trunc, long prec)
{
    long i;
    fmprb_t logc;

    fmprb_init(logc);
    fmprb_log(logc, c, prec);
    fmprb_mul(res + 0, logc, d, prec);
    fmprb_exp(res + 0, res + 0, prec);

    for (i = 1; i < trunc; i++)
    {
        fmprb_mul(res + i, res + i - 1, logc, prec);
        fmprb_div_ui(res + i, res + i, i, prec);
    }

    fmprb_clear(logc);
}

void
_fmprb_poly_zeta_series(fmprb_ptr res, fmprb_srcptr h, long hlen, const fmprb_t a, int deflate, long len, long prec)
{
    long i;
    fmpcb_t cs, ca;
    fmpcb_ptr z;
    fmprb_ptr t, u;

    if (fmprb_contains_nonpositive(a))
    {
        _fmprb_vec_indeterminate(res, len);
        return;
    }

    hlen = FLINT_MIN(hlen, len);

    z = _fmpcb_vec_init(len);
    t = _fmprb_vec_init(len);
    u = _fmprb_vec_init(len);
    fmpcb_init(cs);
    fmpcb_init(ca);

    /* use reflection formula */
    if (fmpr_sgn(fmprb_midref(h)) < 0 && fmprb_is_one(a))
    {
        /* zeta(s) = (2*pi)**s * sin(pi*s/2) / pi * gamma(1-s) * zeta(1-s) */
        fmprb_t pi;
        fmprb_ptr f, s1, s2, s3, s4;

        fmprb_init(pi);
        f = _fmprb_vec_init(2);
        s1 = _fmprb_vec_init(len);
        s2 = _fmprb_vec_init(len);
        s3 = _fmprb_vec_init(len);
        s4 = _fmprb_vec_init(len);

        fmprb_const_pi(pi, prec);

        /* s1 = (2*pi)**s */
        fmprb_mul_2exp_si(pi, pi, 1);
        _fmprb_poly_pow_cpx(s1, pi, h, len, prec);
        fmprb_mul_2exp_si(pi, pi, -1);

        /* s2 = sin(pi*s/2) / pi */
        fmprb_mul_2exp_si(pi, pi, -1);
        fmprb_mul(f, pi, h, prec);
        fmprb_set(f + 1, pi);
        fmprb_mul_2exp_si(pi, pi, 1);
        _fmprb_poly_sin_series(s2, f, 2, len, prec);
        _fmprb_vec_scalar_div(s2, s2, len, pi, prec);

        /* s3 = gamma(1-s) */
        fmprb_sub_ui(f, h, 1, prec);
        fmprb_neg(f, f);
        fmprb_set_si(f + 1, -1);
        _fmprb_poly_gamma_series(s3, f, 2, len, prec);

        /* s4 = zeta(1-s) */
        fmprb_sub_ui(f, h, 1, prec);
        fmprb_neg(f, f);
        fmpcb_set_fmprb(cs, f);
        fmpcb_one(ca);
        zeta_series(z, cs, ca, 0, len, prec);
        for (i = 0; i < len; i++)
            fmprb_set(s4 + i, fmpcb_realref(z + i));
        for (i = 1; i < len; i += 2)
            fmprb_neg(s4 + i, s4 + i);

        _fmprb_poly_mullow(u, s1, len, s2, len, len, prec);
        _fmprb_poly_mullow(s1, s3, len, s4, len, len, prec);
        _fmprb_poly_mullow(t, u, len, s1, len, len, prec);

        /* add 1/(1-(s+t)) = 1/(1-s) + t/(1-s)^2 + ... */
        if (deflate)
        {
            fmprb_sub_ui(u, h, 1, prec);
            fmprb_neg(u, u);
            fmprb_inv(u, u, prec);
            for (i = 1; i < len; i++)
                fmprb_mul(u + i, u + i - 1, u, prec);
            _fmprb_vec_add(t, t, u, len, prec);
        }

        fmprb_clear(pi);
        _fmprb_vec_clear(f, 2);
        _fmprb_vec_clear(s1, len);
        _fmprb_vec_clear(s2, len);
        _fmprb_vec_clear(s3, len);
        _fmprb_vec_clear(s4, len);
    }
    else
    {
        fmpcb_set_fmprb(cs, h);
        fmpcb_set_fmprb(ca, a);
        zeta_series(z, cs, ca, deflate, len, prec);
        for (i = 0; i < len; i++)
            fmprb_set(t + i, fmpcb_realref(z + i));
    }

    /* compose with nonconstant part */
    fmprb_zero(u);
    _fmprb_vec_set(u + 1, h + 1, hlen - 1);
    _fmprb_poly_compose_series(res, t, len, u, hlen, len, prec);

    _fmpcb_vec_clear(z, len);
    _fmprb_vec_clear(t, len);
    _fmprb_vec_clear(u, len);
    fmpcb_init(cs);
    fmpcb_init(ca);
}

void
fmprb_poly_zeta_series(fmprb_poly_t res, const fmprb_poly_t f, const fmprb_t a, int deflate, long n, long prec)
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
        _fmprb_poly_zeta_series(res->coeffs, t, 1, a, deflate, n, prec);
        fmprb_clear(t);
    }
    else
    {
        _fmprb_poly_zeta_series(res->coeffs, f->coeffs, f->length, a, deflate, n, prec);
    }

    _fmprb_poly_set_length(res, n);
    _fmprb_poly_normalise(res);
}

