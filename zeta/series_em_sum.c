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

#include "zeta.h"
#include "fmpcb.h"
#include "fmpcb_poly.h"
#include "bernoulli.h"

/* res = src * (c + x) */
void _fmpcb_poly_mullow_cpx(fmpcb_ptr res, fmpcb_srcptr src, long len, const fmpcb_t c, long trunc, long prec)
{
    long i;

    if (len < trunc)
        fmpcb_set(res + len, src + len - 1);

    for (i = len - 1; i > 0; i--)
    {
        fmpcb_mul(res + i, src + i, c, prec);
        fmpcb_add(res + i, res + i, src + i - 1, prec);
    }

    fmpcb_mul(res, src, c, prec);
}

/* series of 1/N^(c+x) */
void _fmpcb_poly_ui_invpow_cpx(fmpcb_ptr res, ulong N, const fmpcb_t c, long trunc, long prec)
{
    long i;

    fmpcb_t logN;
    fmpcb_init(logN);

    fmprb_log_ui(fmpcb_realref(logN), N, prec);
    fmpcb_mul(res + 0, logN, c, prec);
    fmpcb_neg(res + 0, res + 0);

    fmpcb_exp(res + 0, res + 0, prec);

    for (i = 1; i < trunc; i++)
    {
        fmpcb_mul(res + i, res + i - 1, logN, prec);
        fmpcb_div_si(res + i, res + i, -i, prec);
    }

    fmpcb_clear(logN);
}

void _fmpcb_poly_fmpcb_invpow_cpx(fmpcb_ptr res, const fmpcb_t N, const fmpcb_t c, long trunc, long prec)
{
    long i;
    fmpcb_t logN;

    fmpcb_init(logN);
    fmpcb_log(logN, N, prec);
    fmpcb_mul(res + 0, logN, c, prec);
    fmpcb_neg(res + 0, res + 0);
    fmpcb_exp(res + 0, res + 0, prec);

    for (i = 1; i < trunc; i++)
    {
        fmpcb_mul(res + i, res + i - 1, logN, prec);
        fmpcb_div_si(res + i, res + i, -i, prec);
    }

    fmpcb_clear(logN);
}

static __inline__ int
fmpcb_is_int(const fmpcb_t z)
{
    return fmprb_is_zero(fmpcb_imagref(z)) && fmprb_is_int(fmpcb_realref(z));
}

static __inline__ void
_fmpcb_vec_scalar_div_fmprb(fmpcb_ptr res, fmpcb_srcptr vec, long len, const fmprb_t c, long prec)
{
    long i;
    for (i = 0; i < len; i++)
    {
        fmprb_div(fmpcb_realref(res + i), fmpcb_realref(vec + i), c, prec);
        fmprb_div(fmpcb_imagref(res + i), fmpcb_imagref(vec + i), c, prec);
    }
}

void
zeta_series_em_sum(fmpcb_ptr z, const fmpcb_t s, const fmpcb_t a, int deflate, ulong N, ulong M, long d, long prec)
{
    fmpcb_ptr t, u, v, term, sum;
    fmpcb_t splus, Na, rec;
    fmprb_t x;
    fmpz_t c;
    long i;
    ulong r, n;
    int aint;

    BERNOULLI_ENSURE_CACHED(2 * M);

    t = _fmpcb_vec_init(d + 1);
    u = _fmpcb_vec_init(d);
    v = _fmpcb_vec_init(d);
    term = _fmpcb_vec_init(d);
    sum = _fmpcb_vec_init(d);
    fmpcb_init(splus);
    fmpcb_init(Na);
    fmpcb_init(rec);
    fmprb_init(x);
    fmpz_init(c);

    /* sum 1/(n+a)^(s+x) */
    if (fmpcb_is_one(a) && (d == 1))
    {
        fmpcb_ptr pows;
        long j;

        pows = _fmpcb_vec_init(N + 1);
        fmpcb_one(pows + 1);

        for (i = 2; i <= N; i++)
        {
            if (fmpcb_is_zero(pows + i))
            {
                fmprb_log_ui(fmpcb_realref(pows + i), i, prec);
                fmprb_zero(fmpcb_imagref(pows + i));
                fmpcb_mul(pows + i, pows + i, s, prec);
                fmpcb_neg(pows + i, pows + i);
                fmpcb_exp(pows + i, pows + i, prec);
            }

            for (j = 2; j <= i && i * j <= N; j++)
                if (fmpcb_is_zero(pows + i * j))
                    fmpcb_mul(pows + i * j, pows + i, pows + j, prec);
        }

        for (i = 1; i <= N; i++)
            fmpcb_add(sum, sum, pows + i, prec);

        _fmpcb_vec_clear(pows, N + 1);
    }
    else
    {
        for (n = 0; n < N; n++)
        {
            /* printf("sum 1: %ld %ld\n", n, N); */
            fmpcb_add_ui(Na, a, n, prec);
            _fmpcb_poly_fmpcb_invpow_cpx(t, Na, s, d, prec);
            _fmpcb_vec_add(sum, sum, t, d, prec);
        }
    }

    /* t = 1/(N+a)^(s+x); we might need one extra term for deflation */
    fmpcb_add_ui(Na, a, N, prec);
    _fmpcb_poly_fmpcb_invpow_cpx(t, Na, s, d + 1, prec);

    /* sum += (N+a) * 1/((s+x)-1) * t */
    if (!deflate)
    {
        /* u = 1/(s+x) has series [1/(s-1), -1/(s-1)^2, 1/(s-1)^3, ...] */
        fmpcb_sub_ui(u + 0, s, 1, prec);
        fmpcb_inv(u + 0, u + 0, prec);
        for (i = 1; i < d; i++)
            fmpcb_mul(u + i, u + i - 1, u + 0, prec);
        for (i = 1; i < d; i += 2)
            fmpcb_neg(u + i, u + i);

        _fmpcb_poly_mullow(v, u, d, t, d, d, prec);
        _fmpcb_vec_scalar_mul(v, v, d, Na, prec);
        _fmpcb_vec_add(sum, sum, v, d, prec);
    }
    /* sum += ((N+a)^(1-(s+x)) - 1) / ((s+x) - 1) */
    else
    {
        /* at s = 1, this becomes (N*t - 1)/x, i.e. just remove one coeff  */
        if (fmpcb_is_one(s))
        {
            for (i = 0; i < d; i++)
                fmpcb_mul(u + i, t + i + 1, Na, prec);
            _fmpcb_vec_add(sum, sum, u, d, prec);
        }
        else
        {
            /* TODO: this is numerically unstable for large derivatives,
                and divides by zero if s contains 1. We want a good
                way to evaluate the power series ((N+a)^y - 1) / y where y has
                nonzero constant term, without doing a division.
                How is this best done? */

            _fmpcb_vec_scalar_mul(t, t, d, Na, prec);
            fmpcb_sub_ui(t + 0, t + 0, 1, prec);
            fmpcb_sub_ui(u + 0, s, 1, prec);
            fmpcb_inv(u + 0, u + 0, prec);
            for (i = 1; i < d; i++)
                fmpcb_mul(u + i, u + i - 1, u + 0, prec);
            for (i = 1; i < d; i += 2)
                fmpcb_neg(u + i, u + i);
            _fmpcb_poly_mullow(v, u, d, t, d, d, prec);
            _fmpcb_vec_add(sum, sum, v, d, prec);
            _fmpcb_poly_fmpcb_invpow_cpx(t, Na, s, d, prec);
        }
    }

    /* sum += u = 1/2 * t */
    _fmpcb_vec_scalar_mul_2exp_si(u, t, d, -1L);
    _fmpcb_vec_add(sum, sum, u, d, prec);

    /* term = u * (s+x) / N */
    _fmpcb_poly_mullow_cpx(u, u, d, s, d, prec);
    _fmpcb_vec_scalar_div(term, u, d, Na, prec);

    /* (N+a)^2 or 1/(N+a)^2 */
    fmpcb_mul(Na, Na, Na, prec);
    aint = fmpcb_is_int(Na);

    if (!aint)
        fmpcb_inv(Na, Na, prec);

    for (r = 1; r <= M; r++)
    {
        /* printf("sum 2: %ld %ld\n", r, M); */

        /* sum += bernoulli number * term */
        fmprb_set_round_fmpz(x, fmpq_numref(bernoulli_cache + 2 * r), prec);
        fmprb_div_fmpz(x, x, fmpq_denref(bernoulli_cache + 2 * r), prec);

        _fmpcb_vec_scalar_mul_fmprb(u, term, d, x, prec);
        _fmpcb_vec_add(sum, sum, u, d, prec);

        /* multiply term by ((s+x)+2r-1)((s+x)+2r) / ((N+a)^2 * (2*r+1)*(2*r+2)) */
        fmpcb_set(splus, s);
        fmprb_add_ui(fmpcb_realref(splus), fmpcb_realref(splus), 2*r-1, prec);
        _fmpcb_poly_mullow_cpx(term, term, d, splus, d, prec);
        fmprb_add_ui(fmpcb_realref(splus), fmpcb_realref(splus), 1, prec);
        _fmpcb_poly_mullow_cpx(term, term, d, splus, d, prec);

        /* TODO: combine with previous multiplication? */
        if (aint)
        {
            fmprb_mul_ui(x, fmpcb_realref(Na), 2*r+1, prec);
            fmprb_mul_ui(x, x, 2*r+2, prec);
            _fmpcb_vec_scalar_div_fmprb(term, term, d, x, prec);
        }
        else
        {
            fmpz_set_ui(c, 2*r+1);
            fmpz_mul_ui(c, c, 2*r+2);
            fmpcb_div_fmpz(rec, Na, c, prec);
            _fmpcb_vec_scalar_mul(term, term, d, rec, prec);
        }
    }

    _fmpcb_vec_set(z, sum, d);

    _fmpcb_vec_clear(t, d + 1);
    _fmpcb_vec_clear(u, d);
    _fmpcb_vec_clear(v, d);
    _fmpcb_vec_clear(term, d);
    _fmpcb_vec_clear(sum, d);
    fmpcb_clear(splus);
    fmpcb_clear(Na);
    fmpcb_clear(rec);
    fmprb_clear(x);
    fmpz_clear(c);
}

