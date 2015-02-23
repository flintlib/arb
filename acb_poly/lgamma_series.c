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

#include "acb_poly.h"

void
_acb_log_rising_correct_branch(acb_t t,
        const acb_t t_wrong, const acb_t z, ulong r, long prec);

void acb_gamma_stirling_choose_param(int * reflect, long * r, long * n,
    const acb_t x, int use_reflect, int digamma, long prec);

void
_acb_poly_gamma_stirling_eval(acb_ptr res, const acb_t z, long n, long num, long prec);

static __inline__ void
_log_rising_ui_series(acb_ptr t, const acb_t x, long r, long len, long prec)
{
    acb_struct f[2];
    long rflen;

    acb_init(f);
    acb_init(f + 1);

    acb_set(f, x);
    acb_one(f + 1);

    rflen = FLINT_MIN(len, r + 1);
    _acb_poly_rising_ui_series(t, f, FLINT_MIN(2, len), r, rflen, prec);
    _acb_poly_log_series(t, t, rflen, len, prec);

    _acb_log_rising_correct_branch(t, t, x, r, prec);

    acb_clear(f);
    acb_clear(f + 1);
}

void
_acb_poly_lgamma_series(acb_ptr res, acb_srcptr h, long hlen, long len, long prec)
{
    int reflect;
    long i, r, n, wp;
    acb_t zr;
    acb_ptr t, u;

    hlen = FLINT_MIN(hlen, len);

    if (hlen == 1)
    {
        acb_lgamma(res, h, prec);
        if (acb_is_finite(res))
            _acb_vec_zero(res + 1, len - 1);
        else
            _acb_vec_indeterminate(res + 1, len - 1);
        return;
    }

    if (len == 2)
    {
        acb_t v;
        acb_init(v);
        acb_set(v, h + 1);
        acb_digamma(res + 1, h, prec);
        acb_lgamma(res, h, prec);
        acb_mul(res + 1, res + 1, v, prec);
        acb_clear(v);
        return;
    }

    /* use real code for real input and output */
    if (_acb_vec_is_real(h, hlen) && arb_is_positive(acb_realref(h)))
    {
        arb_ptr tmp = _arb_vec_init(len);
        for (i = 0; i < hlen; i++)
            arb_set(tmp + i, acb_realref(h + i));
        _arb_poly_lgamma_series(tmp, tmp, hlen, len, prec);
        for (i = 0; i < len; i++)
            acb_set_arb(res + i, tmp + i);
        _arb_vec_clear(tmp, len);
        return;
    }

    wp = prec + FLINT_BIT_COUNT(prec);

    t = _acb_vec_init(len);
    u = _acb_vec_init(len);
    acb_init(zr);

    /* use Stirling series */
    acb_gamma_stirling_choose_param(&reflect, &r, &n, h, 0, 0, wp);
    acb_add_ui(zr, h, r, wp);
    _acb_poly_gamma_stirling_eval(u, zr, n, len, wp);

    if (r != 0)
    {
        _log_rising_ui_series(t, h, r, len, wp);
        _acb_vec_sub(u, u, t, len, wp);
    }

    /* compose with nonconstant part */
    acb_zero(t);
    _acb_vec_set(t + 1, h + 1, hlen - 1);
    _acb_poly_compose_series(res, u, len, t, hlen, len, prec);

    acb_clear(zr);
    _acb_vec_clear(t, len);
    _acb_vec_clear(u, len);
}

void
acb_poly_lgamma_series(acb_poly_t res, const acb_poly_t f, long n, long prec)
{
    acb_poly_fit_length(res, n);

    if (f->length == 0 || n == 0)
        _acb_vec_indeterminate(res->coeffs, n);
    else
        _acb_poly_lgamma_series(res->coeffs, f->coeffs, f->length, n, prec);

    _acb_poly_set_length(res, n);
    _acb_poly_normalise(res);
}

