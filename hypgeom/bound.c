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
#include "double_extras.h"
#include "hypgeom.h"

static __inline__ double d_root(double x, int r)
{
    if (r == 1)
        return x;
    if (r == 2)
        return sqrt(x);
    return pow(x, 1. / r);
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
        fmpr_set_round_fmpz(x, t, prec, FMPR_RND_DOWN);
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
        fmpr_set_round_fmpz(x, t, prec, FMPR_RND_UP);
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

/* FIXME: use more than doubles for this */
long hypgeom_root_bound(const fmpr_t z, int r)
{
    if (r == 0)
    {
        return 0;
    }
    else
    {
        double zd;
        zd = fmpr_get_d(z, FMPR_RND_UP);
        return d_root(zd, r) + 2;
    }
}

/*
Given T(K), compute bound for T(n) z^n.

We need to multiply by

z^n * 1/rf(K+1,m)^r * (rf(K+1,m)/rf(K+1-A,m)) * (rf(K+1-B,m)/rf(K+1-2B,m))

where m = n - K. This is equal to

z^n * 

(K+A)! (K-2B)! (K-B+m)!
-----------------------    * ((K+m)! / K!)^(1-r)
(K-B)! (K-A+m)! (K-2B+m)!
*/
void
hypgeom_term_bound(fmpr_t Tn, const fmpr_t TK, long K, long A, long B, int r, const fmpr_t z, long n, long wp)
{
    fmpr_t t, u, num, den;
    long m;

    fmpr_init(t);
    fmpr_init(u);
    fmpr_init(num);
    fmpr_init(den);

    m = n - K;
    if (m < 0)
        abort();

    /* TK * z^n */
    fmpr_pow_sloppy_ui(t, z, n, wp, FMPR_RND_UP);
    fmpr_mul(num, TK, t, wp, FMPR_RND_UP);

    /* (K+A)! (K-2B)! (K-B+m)!, upper bounding */
    fmpr_gamma_ui_ubound(t, K+A+1, wp);
    fmpr_mul(num, num, t, wp, FMPR_RND_UP);
    fmpr_gamma_ui_ubound(t, K-2*B+1, wp);
    fmpr_mul(num, num, t, wp, FMPR_RND_UP);
    fmpr_gamma_ui_ubound(t, K-B+m, wp);
    fmpr_mul(num, num, t, wp, FMPR_RND_UP);

    /* (K-B)! (K-A+m)! (K-2B+m)!, lower bounding */
    fmpr_gamma_ui_lbound(den, K-B+1, wp);
    fmpr_gamma_ui_lbound(t, K-A+m+1, wp);
    fmpr_mul(den, den, t, wp, FMPR_RND_DOWN);
    fmpr_gamma_ui_lbound(t, K-2*B+m+1, wp);
    fmpr_mul(den, den, t, wp, FMPR_RND_DOWN);

    /* ((K+m)! / K!)^(1-r) */
    if (r == 0)
    {
        fmpr_gamma_ui_ubound(t, K+m+1, wp);
        fmpr_mul(num, num, t, wp, FMPR_RND_UP);
        fmpr_gamma_ui_lbound(t, K+1, wp);
        fmpr_mul(den, den, t, wp, FMPR_RND_DOWN);
    }
    else if (r != 1)
    {
        fmpr_gamma_ui_ubound(t, K+1, wp);
        fmpr_gamma_ui_lbound(u, K+m+1, wp);
        fmpr_div(t, t, u, wp, FMPR_RND_UP);
        fmpr_pow_sloppy_ui(t, t, r-1, wp, FMPR_RND_UP);
        fmpr_mul(num, num, t, wp, FMPR_RND_UP);
    }

    fmpr_div(Tn, num, den, wp, FMPR_RND_UP);

    fmpr_clear(t);
    fmpr_clear(u);
    fmpr_clear(num);
    fmpr_clear(den);
}

long
hypgeom_bound(fmpr_t error, int r,
    long A, long B, long K, const fmpr_t TK, const fmpr_t z, long prec)
{
    fmpr_t Tn, t, u, one, tol, num, den;
    long wp = FMPRB_RAD_PREC;
    long n, m;

    fmpr_init(Tn);
    fmpr_init(t);
    fmpr_init(u);
    fmpr_init(one);
    fmpr_init(tol);
    fmpr_init(num);
    fmpr_init(den);

    fmpr_one(one);
    fmpr_set_ui_2exp_si(tol, 1UL, -prec);

    /* approximate number of needed terms */
    n = hypgeom_estimate_terms(z, r, prec);

    /* required for 1 + O(1/k) part to be decreasing */
    n = FLINT_MAX(n, K + 1);

    /* required for z^k / (k!)^r to be decreasing */
    m = hypgeom_root_bound(z, r);
    n = FLINT_MAX(n, m);

    /*  We now have |R(k)| <= G(k) where G(k) is monotonically decreasing,
        and can bound the tail using a geometric series as soon
        as soon as G(k) < 1.  */

    /* bound T(n-1) */
    hypgeom_term_bound(Tn, TK, K, A, B, r, z, n-1, wp);

    while (1)
    {
        /* bound R(n) */
        fmpr_mul_ui(num, z, n, wp, FMPR_RND_UP);
        fmpr_mul_ui(num, num, n - B, wp, FMPR_RND_UP);
        fmpr_set_ui(den, n - A);
        fmpr_mul_ui(den, den, n - 2*B, wp, FMPR_RND_DOWN);

        if (r != 0)
        {
            fmpr_set_ui(u, n);
            fmpr_pow_sloppy_ui(u, u, r, wp, FMPR_RND_DOWN);
            fmpr_mul(den, den, u, wp, FMPR_RND_DOWN);
        }

        fmpr_div(t, num, den, wp, FMPR_RND_UP);

        /* multiply bound for T(n-1) by bound for R(n) to bound T(n) */
        fmpr_mul(Tn, Tn, t, wp, FMPR_RND_UP);

        /* geometric series termination check */
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
        n++;
    }

    fmpr_clear(Tn);
    fmpr_clear(t);
    fmpr_clear(u);
    fmpr_clear(one);
    fmpr_clear(tol);
    fmpr_clear(num);
    fmpr_clear(den);

    return n;
}
