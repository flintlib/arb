/*
    Copyright (C) 2013-2014 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb.h"

static void
bsplit(fmpz_t p1, fmpz_t q1, fmpz_t r1,
        const fmpz_t p, const fmpz_t q, slong a, slong b, int cont)
{
    if (b - a == 1)
    {
        if (a == 0)
        {
            fmpz_set(p1, p);
            fmpz_mul_ui(q1, q, 2*a+1);
            fmpz_mul_ui(r1, p, 2*a+1);
        }
        else
        {
            fmpz_mul(p1, p, p);
            fmpz_mul(q1, q, q);
            fmpz_mul_ui(q1, q1, 2*a+1);
            fmpz_mul_ui(r1, p1, 2*a+1);
        }
    }
    else
    {
        fmpz_t p2, q2, r2;
        slong m;

        m = (a + b) / 2;

        fmpz_init(p2);
        fmpz_init(q2);
        fmpz_init(r2);

        bsplit(p1, q1, r1, p, q, a, m, 1);
        bsplit(p2, q2, r2, p, q, m, b, 1);

        fmpz_mul(p1, p1, q2);
        fmpz_addmul(p1, r1, p2);
        fmpz_mul(q1, q1, q2);
        if (cont)
            fmpz_mul(r1, r1, r2);

        fmpz_clear(p2);
        fmpz_clear(q2);
        fmpz_clear(r2);

    }
}

#define LOG2 0.69314718055994530942

/* Assumption: p/q <= 2 */
static void
atanh_bs(arb_t s, ulong p, ulong q, slong prec)
{
    fmpz_t pp, qq, P, Q, R;
    double logqp;
    slong N;

    if (p == 0)
    {
        arb_zero(s);
        return;
    }

    fmpz_init(pp);
    fmpz_init(qq);
    fmpz_init(P);
    fmpz_init(Q);
    fmpz_init(R);

    /* If N >= 1 and p/q <= 1/2, the error is bounded by (p/q)^(2N+1).
    For error <= 2^-prec, it is sufficient to pick
    N >= (1/2) * (prec * log(2) / log(q/p) - 1). */
    logqp = mag_d_log_lower_bound(q / (double) p) * (1.0 - 1e-12);
    N = ceil((prec * (0.5 * LOG2) / logqp) * (1.0 + 1e-12));
    N = FLINT_MAX(N, 1);

    fmpz_set_ui(pp, p);
    fmpz_set_ui(qq, q);

    bsplit(P, Q, R, pp, qq, 0, N, 0);

    arb_fmpz_div_fmpz(s, P, Q, prec);
    arb_add_error_2exp_si(s, -prec);

    fmpz_clear(pp);
    fmpz_clear(qq);
    fmpz_clear(P);
    fmpz_clear(Q);
    fmpz_clear(R);
}

static int n_width(ulong k)
{
    int a, b;
    count_leading_zeros(a, k);
    count_trailing_zeros(b, k);
    return FLINT_BITS - a - b;
}

void
arb_log_ui_from_prev(arb_t s, ulong k, arb_t log_prev, ulong prev, slong prec)
{
    if (prev < 2 || prec < 600 ||
        (prec < ARB_LOG_TAB2_PREC - 64 && n_width(k) <= ARB_LOG_TAB21_BITS + 1)
        || k < prev || (k + prev) < prev ||
        (k - prev) >= 0.25 * (k + prev))
    {
        arb_log_ui(s, k, prec);
    }
    else
    {
        arb_t t;
        ulong p, q;

        arb_init(t);

        p = k - prev;
        q = k + prev;

        if ((p % 2 == 0) && (q % 2 == 0))
        {
            p >>= 1;
            q >>= 1;
        }

        atanh_bs(t, p, q, prec);
        arb_mul_2exp_si(t, t, 1);
        arb_add(s, log_prev, t, prec);

        arb_clear(t);
    }
}

