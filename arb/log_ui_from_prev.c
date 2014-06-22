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

    Copyright (C) 2013-2014 Fredrik Johansson

******************************************************************************/

#include "arb.h"

static void
bsplit(fmpz_t p1, fmpz_t q1, fmpz_t r1,
        const fmpz_t p, const fmpz_t q, long a, long b, int cont)
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
        long m;

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
atanh_bs(arb_t s, ulong p, ulong q, long prec)
{
    fmpz_t pp, qq, P, Q, R;
    double logqp;
    long N;

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

/* Assumption: p/q <= 2 */
static void
atanh_rec(arb_t z, ulong p, ulong q, long prec)
{
    fmpz_t u, t, s, p2, q2, r;
    ulong k;
    mp_limb_t err;

    if (p == 0)
    {
        arb_zero(z);
        return;
    }

    fmpz_init(u);
    fmpz_init(t);
    fmpz_init(s);
    fmpz_init(p2);
    fmpz_init(q2);
    fmpz_init(r);

    prec += 4 + FLINT_BIT_COUNT(prec);

    /* precompute p^2, q^2 */
    fmpz_set_ui(p2, p);
    fmpz_mul_ui(p2, p2, p);
    fmpz_set_ui(q2, q);
    fmpz_mul_ui(q2, q2, q);

    /* s = t = p / q */
    fmpz_set_ui(t, p);
    fmpz_mul_2exp(t, t, prec);
    fmpz_tdiv_q_ui(t, t, q);
    fmpz_set(s, t);

    /* err(s) = 1 ulp */
    err = 1;

    k = 3;
    while (!fmpz_is_zero(t))
    {
        /* t = t * p^2 / q^2 */
        /* invariant: err(t) <= 2 */
        fmpz_mul(t, t, p2);
        fmpz_tdiv_q(t, t, q2);

        /* err(u) <= 2 */
        fmpz_tdiv_q_ui(u, t, k);

        /* s = s + t, plus at most 2 ulp error */
        fmpz_add(s, s, u);
        err += 2;

        k += 2;
    }

    /* The truncation error is less than twice
       the first omitted term, i.e. less than 2*2 ulp */
    err += 4;

    arb_set_fmpz(z, s);
    mag_set_ui(arb_radref(z), err);
    arb_mul_2exp_si(z, z, -prec);

    fmpz_clear(u);
    fmpz_clear(t);
    fmpz_clear(s);
    fmpz_clear(p2);
    fmpz_clear(q2);
    fmpz_clear(r);
}

void
arb_log_ui_from_prev(arb_t s, ulong k, arb_t log_prev, ulong prev, long prec)
{
    if (prev >= 2 &&
        (k + prev) > prev &&               /* no overflow */
        (k - prev) < 0.25 * (k + prev) &&  /* rapid convergence */
        k >= prev)
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

        if (prec <= 256)
            atanh_rec(t, p, q, prec);
        else
            atanh_bs(t, p, q, prec);

        arb_mul_2exp_si(t, t, 1);
        arb_add(s, log_prev, t, prec);

        arb_clear(t);
    }
    else
    {
        arb_log_ui(s, k, prec);
    }
}

