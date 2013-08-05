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

#include "zeta.h"
#include "fmpcb.h"
#include "fmpcb_poly.h"

void _fmpcb_poly_fmpcb_invpow_cpx(fmpcb_ptr res,
    const fmpcb_t N, const fmpcb_t c, long trunc, long prec);

#define POWER(_k) (powers + (((_k)-1)/2) * (len))
#define DIVISOR(_k) (divisors[((_k)-1)/2])

#define COMPUTE_POWER(t, k, kprev) \
  do { \
    if (integer) \
    { \
        fmprb_neg(w, fmpcb_realref(s)); \
        fmprb_set_ui(v, k); \
        fmprb_pow(fmpcb_realref(t), v, w, prec); \
        fmprb_zero(fmpcb_imagref(t)); \
        if (len != 1) \
        { \
            zeta_log_ui_from_prev(logk, k, logk, kprev, prec); \
            kprev = k; \
            fmprb_neg(logk, logk); \
        } \
    } \
    else \
    { \
        zeta_log_ui_from_prev(logk, k, logk, kprev, prec); \
        kprev = k; \
        fmprb_neg(logk, logk); \
        fmprb_mul(w, logk, fmpcb_imagref(s), prec); \
        fmprb_sin_cos(fmpcb_imagref(t), fmpcb_realref(t), w, prec); \
        if (critical_line) \
        { \
            fmprb_rsqrt_ui(w, k, prec); \
            fmpcb_mul_fmprb(t, t, w, prec); \
        } \
        else \
        { \
            fmprb_mul(w, fmpcb_realref(s), logk, prec); \
            fmprb_exp(w, w, prec); \
            fmpcb_mul_fmprb(t, t, w, prec); \
        } \
    } \
    for (i = 1; i < len; i++) \
    { \
        fmpcb_mul_fmprb(t + i, t + i - 1, logk, prec); \
        fmpcb_div_ui(t + i, t + i, i, prec); \
    } \
    fmprb_neg(logk, logk); \
  } while (0); \

void
zeta_powsum_one_series_sieved(fmpcb_ptr z, const fmpcb_t s, long n, long len, long prec)
{
    long * divisors;
    long powers_alloc;
    long i, j, k, kprev, power_of_two, horner_point;
    int critical_line, integer;

    fmpcb_ptr powers;
    fmpcb_ptr t, u, x;
    fmpcb_ptr p1, p2;
    fmprb_t logk, v, w;

    critical_line = fmprb_is_exact(fmpcb_realref(s)) &&
        (fmpr_cmp_2exp_si(fmprb_midref(fmpcb_realref(s)), -1) == 0);

    integer = fmprb_is_zero(fmpcb_imagref(s)) && fmprb_is_int(fmpcb_realref(s));

    divisors = flint_calloc(n / 2 + 1, sizeof(long));
    powers_alloc = (n / 6 + 1) * len;
    powers = _fmpcb_vec_init(powers_alloc);

    for (i = 3; i <= n; i += 2)
        for (j = 3 * i; j <= n; j += 2 * i)
            DIVISOR(j) = i;

    t = _fmpcb_vec_init(len);
    u = _fmpcb_vec_init(len);
    x = _fmpcb_vec_init(len);
    fmprb_init(logk);
    fmprb_init(v);
    fmprb_init(w);

    power_of_two = 1;
    while (power_of_two * 2 <= n)
        power_of_two *= 2;
    horner_point = n / power_of_two;

    _fmpcb_vec_zero(z, len);

    kprev = 0;
    COMPUTE_POWER(x, 2, kprev);

    for (k = 1; k <= n; k += 2)
    {
        /* t = k^(-s) */
        if (DIVISOR(k) == 0)
        {
            COMPUTE_POWER(t, k, kprev);
        }
        else
        {
            p1 = POWER(DIVISOR(k));
            p2 = POWER(k / DIVISOR(k));

            if (len == 1)
                fmpcb_mul(t, p1, p2, prec);
            else
                _fmpcb_poly_mullow(t, p1, len, p2, len, len, prec);
        }

        if (k * 3 <= n)
            _fmpcb_vec_set(POWER(k), t, len);

        _fmpcb_vec_add(u, u, t, len, prec);

        while (k == horner_point && power_of_two != 1)
        {
            _fmpcb_poly_mullow(t, z, len, x, len, len, prec);
            _fmpcb_vec_add(z, t, u, len, prec);

            power_of_two /= 2;
            horner_point = n / power_of_two;
            horner_point -= (horner_point % 2 == 0);
        }
    }

    _fmpcb_poly_mullow(t, z, len, x, len, len, prec);
    _fmpcb_vec_add(z, t, u, len, prec);

    flint_free(divisors);
    _fmpcb_vec_clear(powers, powers_alloc);
    _fmpcb_vec_clear(t, len);
    _fmpcb_vec_clear(u, len);
    _fmpcb_vec_clear(x, len);
    fmprb_clear(logk);
    fmprb_clear(v);
    fmprb_clear(w);
}

