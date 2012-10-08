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

#include <math.h>
#include "fmprb.h"
#include "double_extras.h"

double fmpr_get_d(const fmpr_t x);

#define LOG2 0.69314718055994530942
#define EXP1 2.7182818284590452354

double d_root(double x, int r)
{
    if (r == 1)
        return x;
    if (r == 2)
        return sqrt(x);
    return pow(x, 1. / r);
}

/*
Estimate the truncation point to obtain accuracy 2^(-prec) with the
hypergeometric series |z|^k / (k!)^r.
*/
long
estimate_nterms(double z, int r, long prec)
{
    double y;

    z = fabs(z);

    if (z == 0)
        return 1;

    if (r == 0)
    {
        if (z >= 1)
        {
            printf("z must be smaller than 1\n");
            abort();
        }

        y = (log(1-z) - prec * LOG2) / log(z) + 1;
    }
    else
    {
        /* Solve k*log(z) - r*(k*log(k)-k) = -prec*log(2)  */
        y = prec * LOG2 / (d_root(z, r) * EXP1 * r);
        y =  prec * LOG2 / (r * d_lambertw(y)) + 1;
    }

    if (y >= LONG_MAX / 2)
    {
        printf("error: series will converge too slowly\n");
        abort();
    }

    return y;
}

void
fmpr_gamma_ui_lbound(fmpr_t x, ulong n, long prec)
{
    if (n == 0) abort();

    if (n < 250)
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_fac_ui(t, n - 1);
        fmpr_set_fmpz(x, t);
        fmpr_set_round(x, x, prec, FMPR_RND_DOWN);
        fmpz_clear(t);
    }
    else
    {
        /* (2 pi/x)^(1/2) * (x/e)^x < Gamma(x) */
        fmpr_t t, u;

        fmpr_init(t);
        fmpr_init(u);

        /* lower bound for 2 pi */
        fmpr_set_ui_2exp_si(t, 843314855, -27);
        fmpr_div_ui(t, t, n, prec, FMPR_RND_DOWN);
        fmpr_sqrt(t, t, prec, FMPR_RND_DOWN);

        /* lower bound for 1/e */
        fmpr_set_ui_2exp_si(u, 197503771, -29);
        fmpr_mul_ui(u, u, n, prec, FMPR_RND_DOWN);

        fmpr_pow_sloppy_ui(u, u, n, prec, FMPR_RND_DOWN);

        fmpr_mul(x, t, u, prec, FMPR_RND_DOWN);

        fmpr_clear(t);
        fmpr_clear(u);
    }
}

void
fmpr_gamma_ui_ubound(fmpr_t x, ulong n, long prec)
{
    if (n == 0) abort();

    if (n < 250)
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_fac_ui(t, n - 1);
        fmpr_set_fmpz(x, t);
        fmpr_set_round(x, x, prec, FMPR_RND_UP);
        fmpz_clear(t);
    }
    else
    {
        fmpr_t t, u;
        fmpr_init(t);

        /* Gamma(x) < e * (x / e)^x -- TODO: use a tighter bound */

        fmpr_init(t);
        fmpr_init(u);

        /* upper bound for 1/e */
        fmpr_set_ui_2exp_si(u, 197503773, -29);
        fmpr_mul_ui(u, u, n, prec, FMPR_RND_UP);
        fmpr_pow_sloppy_ui(u, u, n, prec, FMPR_RND_UP);

        /* upper bound for e */
        fmpr_set_ui_2exp_si(t, 364841613, -27);
        fmpr_mul(x, t, u, prec, FMPR_RND_UP);

        fmpr_clear(t);
        fmpr_clear(u);
    }
}

/* FIXME: assumes no overflow when computing n + r */
static void
fmpr_rfac_uiui_ubound(fmpr_t x, ulong n, ulong r, long prec)
{
    if (r == 0)
    {
        fmpr_one(x);
    }
    else if (r == 1)
    {
        fmpr_set_ui(x, n);
    }
    else
    {
        fmpr_t t;
        fmpr_init(t);
        fmpr_gamma_ui_ubound(x, n + r, prec);
        fmpr_gamma_ui_lbound(t, n, prec);
        fmpr_div(x, x, t, prec, FMPR_RND_UP);
        fmpr_clear(t);
    }
}

static void
fmpr_rfac_uiui_lbound(fmpr_t x, ulong n, ulong r, long prec)
{
    if (r == 0)
    {
        fmpr_one(x);
    }
    else if (r == 1)
    {
        fmpr_set_ui(x, n);
    }
    else
    {
        fmpr_t t;
        fmpr_init(t);
        fmpr_gamma_ui_lbound(x, n + r, prec);
        fmpr_gamma_ui_ubound(t, n, prec);
        fmpr_div(x, x, t, prec, FMPR_RND_DOWN);
        fmpr_clear(t);
    }
}


/*

The general term T(k) is z^k / (k!)^r * R(k) where R(k) = 1 + O(1/k).

We have precomputed integers A, B, K such that for all k > K,
|R(k)| <= F(k) = k(k-B)/((k-A)(k-2B)) = (1+A/(k-A))(1+B/(k-2B)).

*/
long
hypgeom_bound(fmpr_t error, int r,
    long A, long B, long K, const fmpr_t TK, const fmpr_t z, long prec)
{
    fmpr_t Tn, t, u, one, tol;
    long wp = FMPRB_RAD_PREC;
    long n;
    double zd;

    fmpr_init(Tn);
    fmpr_init(t);
    fmpr_init(u);
    fmpr_init(one);
    fmpr_init(tol);

    fmpr_one(one);
    fmpr_set_ui_2exp_si(tol, 1UL, -prec);

    zd = fmpr_get_d(z);

    n = estimate_nterms(zd, r, prec);

    /* required for 1 + O(1/k) part to be decreasing */
    n = FLINT_MAX(n, K + 1);

    /* required for z^k / (k!)^r to be decreasing
       (TODO: don't use doubles for this) */
    if (r > 0)
    {
        long nbd = d_root(zd, r) + 2;
        n = FLINT_MAX(n, nbd);
    }

    /* We are now sure that |R(k)| is either decreasing or strictly
       smaller than 1 for k >= n, which means that we can bound the tail
       using a geometric series as soon as soon as |R(k)| < 1  */

    /* Compute an upper bound for T(n) */
    /* (z^n) / (n!)^r * TK * [(K+1)(K+2)...(n)] * [(K-B+1)(K-B+2)...(n-B)]
                             ---------------------------------------------
                             [(K-A+1)(K-A+2)...(n)] * [(K-2B+1)(K-2B+2)...(n-2B)]
    */

    /* z^n * TK */
    fmpr_pow_sloppy_ui(Tn, z, n, wp, FMPR_RND_UP);
    fmpr_mul(Tn, Tn, TK, wp, FMPR_RND_UP);

    /* divide by (n!)^r */
    if (r != 0)
    {
        fmpr_gamma_ui_lbound(t, n + 1, wp);
        fmpr_ui_div(t, 1UL, t, wp, FMPR_RND_UP);
        fmpr_pow_sloppy_ui(t, t, r, wp, FMPR_RND_UP);
        fmpr_mul(Tn, Tn, t, wp, FMPR_RND_UP);
    }

    fmpr_rfac_uiui_ubound(t, K+1, n-K, wp);
    fmpr_mul(Tn, Tn, t, wp, FMPR_RND_UP);

    fmpr_rfac_uiui_ubound(t, K-B+1, n-K, wp);
    fmpr_mul(Tn, Tn, t, wp, FMPR_RND_UP);

    fmpr_rfac_uiui_lbound(t, K-A+1, n-K, wp);
    fmpr_div(Tn, Tn, t, wp, FMPR_RND_UP);

    fmpr_rfac_uiui_lbound(t, K-2*B+1, n-K, wp);
    fmpr_div(Tn, Tn, t, wp, FMPR_RND_UP);

    while (1)
    {
        /* bound for term ratio: z * F(n) / n^r */

        /* F(n) <= n (n-B) / ((n-A) (n-2B)) */
        fmpr_set_ui(t, n);
        fmpr_mul_ui(t, t, n - B, wp, FMPR_RND_UP);
        fmpr_div_ui(t, t, n - A, wp, FMPR_RND_UP);
        fmpr_div_ui(t, t, n - 2*B, wp, FMPR_RND_UP);

        fmpr_mul(t, t, z, wp, FMPR_RND_UP);

        if (r != 0)
        {
            fmpr_div_ui(u, one, n, wp, FMPR_RND_UP);
            fmpr_pow_sloppy_ui(u, u, r, wp, FMPR_RND_UP);
            fmpr_mul(t, t, u, wp, FMPR_RND_UP);
        }

        /* bound by geometric series: Tn / (1 - t) */
        /* where the term ratio must be < 1 */
        fmpr_sub(u, one, t, wp, FMPR_RND_DOWN);

        if (fmpr_sgn(u) > 0)
        {
            fmpr_div(u, Tn, u, wp, FMPR_RND_UP);

            if (fmpr_cmp(u, tol) < 0)
            {
                fmpr_set(error, u);
                break;
            }
        }

        /* move on to next term */
        fmpr_mul(Tn, Tn, t, wp, FMPR_RND_UP);
        n++;
    }

    fmpr_clear(Tn);
    fmpr_clear(t);
    fmpr_clear(u);
    fmpr_clear(one);
    fmpr_clear(tol);

    return n;
}
