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

    Copyright (C) 2012-2014 Fredrik Johansson

******************************************************************************/

#include "acb_poly.h"

#define POWER(_k) (powers + (((_k)-1)/2) * (len))
#define DIVISOR(_k) (divisors[((_k)-1)/2])

#define COMPUTE_POWER(t, k, kprev) \
  do { \
    if (integer) \
    { \
        arb_neg(w, acb_realref(s)); \
        arb_set_ui(v, k); \
        arb_pow(acb_realref(t), v, w, prec); \
        arb_zero(acb_imagref(t)); \
        if (len != 1) \
        { \
            arb_log_ui_from_prev(logk, k, logk, kprev, prec); \
            kprev = k; \
            arb_neg(logk, logk); \
        } \
    } \
    else \
    { \
        arb_log_ui_from_prev(logk, k, logk, kprev, prec); \
        kprev = k; \
        arb_neg(logk, logk); \
        arb_mul(w, logk, acb_imagref(s), prec); \
        arb_sin_cos(acb_imagref(t), acb_realref(t), w, prec); \
        if (critical_line) \
        { \
            arb_rsqrt_ui(w, k, prec); \
            acb_mul_arb(t, t, w, prec); \
        } \
        else \
        { \
            arb_mul(w, acb_realref(s), logk, prec); \
            arb_exp(w, w, prec); \
            acb_mul_arb(t, t, w, prec); \
        } \
    } \
    for (i = 1; i < len; i++) \
    { \
        acb_mul_arb(t + i, t + i - 1, logk, prec); \
        acb_div_ui(t + i, t + i, i, prec); \
    } \
    arb_neg(logk, logk); \
  } while (0); \

void
_acb_poly_powsum_one_series_sieved(acb_ptr z, const acb_t s, long n, long len, long prec)
{
    long * divisors;
    long powers_alloc;
    long i, j, k, ibound, kprev, power_of_two, horner_point;
    int critical_line, integer;

    acb_ptr powers;
    acb_ptr t, u, x;
    acb_ptr p1, p2;
    arb_t logk, v, w;

    critical_line = arb_is_exact(acb_realref(s)) &&
        (arf_cmp_2exp_si(arb_midref(acb_realref(s)), -1) == 0);

    integer = arb_is_zero(acb_imagref(s)) && arb_is_int(acb_realref(s));

    divisors = flint_calloc(n / 2 + 1, sizeof(long));
    powers_alloc = (n / 6 + 1) * len;
    powers = _acb_vec_init(powers_alloc);

    ibound = n_sqrt(n);
    for (i = 3; i <= ibound; i += 2)
        if (DIVISOR(i) == 0)
            for (j = i * i; j <= n; j += 2 * i)
                DIVISOR(j) = i;

    t = _acb_vec_init(len);
    u = _acb_vec_init(len);
    x = _acb_vec_init(len);
    arb_init(logk);
    arb_init(v);
    arb_init(w);

    power_of_two = 1;
    while (power_of_two * 2 <= n)
        power_of_two *= 2;
    horner_point = n / power_of_two;

    _acb_vec_zero(z, len);

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
                acb_mul(t, p1, p2, prec);
            else
                _acb_poly_mullow(t, p1, len, p2, len, len, prec);
        }

        if (k * 3 <= n)
            _acb_vec_set(POWER(k), t, len);

        _acb_vec_add(u, u, t, len, prec);

        while (k == horner_point && power_of_two != 1)
        {
            _acb_poly_mullow(t, z, len, x, len, len, prec);
            _acb_vec_add(z, t, u, len, prec);

            power_of_two /= 2;
            horner_point = n / power_of_two;
            horner_point -= (horner_point % 2 == 0);
        }
    }

    _acb_poly_mullow(t, z, len, x, len, len, prec);
    _acb_vec_add(z, t, u, len, prec);

    flint_free(divisors);
    _acb_vec_clear(powers, powers_alloc);
    _acb_vec_clear(t, len);
    _acb_vec_clear(u, len);
    _acb_vec_clear(x, len);
    arb_clear(logk);
    arb_clear(v);
    arb_clear(w);
}

