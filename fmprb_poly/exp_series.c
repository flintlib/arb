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

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "fmprb_poly.h"

#define NEWTON_EXP_CUTOFF 60


/* with inverse=1 simultaneously computes g = exp(-x) to length n
with inverse=0 uses g as scratch space, computing
g = exp(-x) only to length (n+1)/2 */
static void
_fmprb_poly_exp_series_newton(fmprb_ptr f, fmprb_ptr g,
    fmprb_srcptr h, long len, long prec, int inverse, long cutoff)
{
    long alloc;
    fmprb_ptr T, U, hprime;

    alloc = 3 * len;
    T = _fmprb_vec_init(alloc);
    U = T + len;
    hprime = U + len;

    _fmprb_poly_derivative(hprime, h, len, prec);
    fmprb_zero(hprime + len - 1);

    NEWTON_INIT(cutoff, len)

    /* f := exp(h) + O(x^m), g := exp(-h) + O(x^m2) */
    NEWTON_BASECASE(n)
    _fmprb_poly_exp_series_basecase(f, h, n, n, prec);
    _fmprb_poly_inv_series(g, f, (n + 1) / 2, (n + 1) / 2, prec);
    NEWTON_END_BASECASE

    /* extend from length m to length n */
    NEWTON_LOOP(m, n)

    long m2 = (m + 1) / 2;
    long l = m - 1; /* shifted for derivative */

    /* g := exp(-h) + O(x^m) */
    _fmprb_poly_mullow(T, f, m, g, m2, m, prec);
    _fmprb_poly_mullow(g + m2, g, m2, T + m2, m - m2, m - m2, prec);
    _fmprb_vec_neg(g + m2, g + m2, m - m2);

    /* U := h' + g (f' - f h') + O(x^(n-1))
        Note: should replace h' by h' mod x^(m-1) */
    _fmprb_vec_zero(f + m, n - m);
    _fmprb_poly_mullow(T, f, n, hprime, n, n, prec); /* should be mulmid */
    _fmprb_poly_derivative(U, f, n, prec); fmprb_zero(U + n - 1); /* should skip low terms */
    _fmprb_vec_sub(U + l, U + l, T + l, n - l, prec);
    _fmprb_poly_mullow(T + l, g, n - m, U + l, n - m, n - m, prec);
    _fmprb_vec_add(U + l, hprime + l, T + l, n - m, prec);

    /* f := f + f * (h - int U) + O(x^n) = exp(h) + O(x^n) */
    _fmprb_poly_integral(U, U, n, prec); /* should skip low terms */
    _fmprb_vec_sub(U + m, h + m, U + m, n - m, prec);
    _fmprb_poly_mullow(f + m, f, n - m, U + m, n - m, n - m, prec);

    /* g := exp(-h) + O(x^n) */
    /* not needed if we only want exp(x) */
    if (n == len && inverse)
    {
        _fmprb_poly_mullow(T, f, n, g, m, n, prec);
        _fmprb_poly_mullow(g + m, g, m, T + m, n - m, n - m, prec);
        _fmprb_vec_neg(g + m, g + m, n - m);
    }

    NEWTON_END_LOOP

    NEWTON_END

    _fmprb_vec_clear(T, alloc);
}

void
_fmprb_poly_exp_series(fmprb_ptr f, fmprb_srcptr h, long hlen, long n, long prec)
{
    hlen = FLINT_MIN(hlen, n);

    if (hlen == 1)
    {
        fmprb_exp(f, h, prec);
        _fmprb_vec_zero(f + 1, n - 1);
    }
    else if (n == 2)
    {
        fmprb_exp(f, h, prec);
        fmprb_mul(f + 1, f, h + 1, prec);  /* safe since hlen >= 2 */
        _fmprb_vec_zero(f + 2, n - 2);
    }
    else if (_fmprb_vec_is_zero(h + 1, hlen - 2)) /* h = a + bx^c */
    {
        long i, j, d = hlen - 1;
        fmprb_t t;
        fmprb_init(t);
        fmprb_set(t, h + d);
        fmprb_exp(f, h, prec);
        for (i = 1, j = d; j < n; j += d, i++)
        {
            fmprb_mul(f + j, f + j - d, t, prec);
            fmprb_div_ui(f + j, f + j, i, prec);
            _fmprb_vec_zero(f + j - d + 1, hlen - 2);
        }
        _fmprb_vec_zero(f + j - d + 1, n - (j - d + 1));
        fmprb_clear(t);
    }
    else if (hlen < NEWTON_EXP_CUTOFF)
    {
        _fmprb_poly_exp_series_basecase(f, h, hlen, n, prec);
    }
    else
    {
        fmprb_ptr g, t;
        fmprb_t u;
        int fix;

        g = _fmprb_vec_init((n + 1) / 2);
        fix = (hlen < n || h == f || !fmprb_is_zero(h));

        if (fix)
        {
            t = _fmprb_vec_init(n);
            _fmprb_vec_set(t + 1, h + 1, hlen - 1);
        }
        else
            t = (fmprb_ptr) h;

        fmprb_init(u);
        fmprb_exp(u, h, prec);

        _fmprb_poly_exp_series_newton(f, g, t, n, prec, 0, NEWTON_EXP_CUTOFF);

        if (!fmprb_is_one(u))
            _fmprb_vec_scalar_mul(f, f, n, u, prec);

        _fmprb_vec_clear(g, (n + 1) / 2);
        if (fix)
            _fmprb_vec_clear(t, n);
        fmprb_clear(u);
    }
}

void
fmprb_poly_exp_series(fmprb_poly_t f, const fmprb_poly_t h, long n, long prec)
{
    long hlen = h->length;

    if (n == 0)
    {
        fmprb_poly_zero(f);
        return;
    }

    if (hlen == 0)
    {
        fmprb_poly_one(f);
        return;
    }

    if (hlen == 1)
        n = 1;

    fmprb_poly_fit_length(f, n);
    _fmprb_poly_exp_series(f->coeffs, h->coeffs, hlen, n, prec);
    _fmprb_poly_set_length(f, n);
    _fmprb_poly_normalise(f);
}
