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

#include "gamma.h"
#include "fmpcb_poly.h"

static void
bsplit(fmpcb_ptr Q, fmpcb_ptr T, const fmpcb_t z, long a, long b, long num, long prec)
{
    if (b - a == 1)
    {
        gamma_stirling_coeff(fmpcb_realref(T), a, 0, prec);
        fmprb_zero(fmpcb_imagref(T));

        if (a == 1)
        {   /* (z + t) */
            fmpcb_set(Q, z);
            if (num > 1) fmpcb_one(Q + 1);
            if (num > 2) fmpcb_zero(Q + 2);
        }             
        else
        {   /* (z + t)^2 */
            fmpcb_mul(Q, z, z, prec);  /* TODO: precompute */
            if (num > 1) fmpcb_mul_2exp_si(Q + 1, z, 1);
            if (num > 2) fmpcb_one(Q + 2);
        }
    }
    else
    {
        long m, n1, n2, q1len, q2len, t1len, t2len, qlen, tlen, alloc;
        fmpcb_ptr Q1, T1, Q2, T2;

        m = a + (b - a) / 2;

        n1 = m - a;
        n2 = b - m;
        q1len = FLINT_MIN(2 * n1 + 1, num);
        t1len = FLINT_MIN(2 * n1 - 1, num);
        q2len = FLINT_MIN(2 * n2 + 1, num);
        t2len = FLINT_MIN(2 * n2 - 1, num);
        qlen = FLINT_MIN(q1len + q2len - 1, num);
        tlen = FLINT_MIN(t1len + q2len - 1, num);

        alloc = q1len + q2len + t1len + t2len;
        Q1 = _fmpcb_vec_init(alloc);
        Q2 = Q1 + q1len;
        T1 = Q2 + q2len;
        T2 = T1 + t1len;

        bsplit(Q1, T1, z, a, m, num, prec);
        bsplit(Q2, T2, z, m, b, num, prec);

        _fmpcb_poly_mullow(Q, Q2, q2len, Q1, q1len, qlen, prec);
        _fmpcb_poly_mullow(T, Q2, q2len, T1, t1len, tlen, prec);
        _fmpcb_poly_add(T, T, tlen, T2, t2len, prec);

        _fmpcb_vec_clear(Q1, alloc);
    }
}

void
_fmpcb_poly_mullow_cpx(fmpcb_ptr res, fmpcb_srcptr src, long len, const fmpcb_t c, long trunc, long prec);

void
_fmpcb_poly_log_cpx_series(fmpcb_ptr res, const fmpcb_t c, long num, long prec)
{
    long i;

    for (i = 0; i < num; i++)
    {
        if (i == 0)
            fmpcb_log(res + i, c, prec);
        else if (i == 1)
            fmpcb_inv(res + i, c, prec);
        else
            fmpcb_mul(res + i, res + i - 1, res + 1, prec);
    }

    for (i = 2; i < num; i++)
    {
        fmpcb_div_ui(res + i, res + i, i, prec);

        if (i % 2 == 0)
            fmpcb_neg(res + i, res + i);
    }
}

void
gamma_stirling_eval_fmpcb_series(fmpcb_ptr res, const fmpcb_t z, long n, long num, long prec)
{
    long tlen, qlen;
    fmpcb_ptr T, Q;
    fmpr_ptr err;
    fmpcb_t c;

    T = _fmpcb_vec_init(num);
    Q = _fmpcb_vec_init(num);
    err = _fmpr_vec_init(num);
    fmpcb_init(c);

    gamma_stirling_bound_fmpcb(err, z, 0, num, n);

    if (n <= 1)
    {
        _fmpcb_vec_zero(res, num);
    }
    else
    {
        qlen = FLINT_MIN(2 * (n - 1) + 1, num);
        tlen = FLINT_MIN(2 * (n - 1) - 1, num);
        bsplit(Q, T, z, 1, n, num, prec);
        _fmpcb_poly_div_series(res, T, tlen, Q, qlen, num, prec);
    }

    /* ((z-1/2) + t) * log(z+t) */
    _fmpcb_poly_log_cpx_series(T, z, num, prec);
    fmpcb_one(c);
    fmpcb_mul_2exp_si(c, c, -1);
    fmpcb_sub(c, z, c, prec);
    _fmpcb_poly_mullow_cpx(T, T, num, c, num, prec);

    /* constant term */
    fmprb_const_log_sqrt2pi(fmpcb_realref(c), prec);
    fmprb_zero(fmpcb_imagref(c));
    fmpcb_add(T, T, c, prec);

    /* subtract (z+t) */
    fmpcb_sub(T, T, z, prec);
    if (num > 1)
        fmpcb_sub_ui(T + 1, T + 1, 1, prec);

    _fmpcb_vec_add(res, res, T, num, prec);

    _fmpcb_vec_add_error_fmpr_vec(res, err, num);

    _fmpcb_vec_clear(T, num);
    _fmpcb_vec_clear(Q, num);
    _fmpr_vec_clear(err, num);
    fmpcb_clear(c);
}

