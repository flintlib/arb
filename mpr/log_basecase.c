/*=============================================================================

    This file is part of ARB.

    ARB is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    ARB is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with ARB; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include <alloca.h>
#include "mpr.h"

/* assumed accuracy of initial approximation
   (the current value is quite conservative) */
#define INITIAL_PREC 48

mp_limb_t
mpr_log_initial(mp_limb_t x)
{
    double t;

    t = ((double) x) * 5.42101086242752217e-20;   /* 2^(-64) */
    t = log(1.0 + t);
    return (mp_limb_t) (t * 18446744073709551616.0);  /* 2^64 */
}

void
mpr_log_basecase(mp_ptr y, mp_srcptr x, long limbs, mp_bitcnt_t tol_bits)
{
    /* TODO: dynamic allocation */
    mp_limb_t t[16], u[16], w[16];
    long n, i, bit_steps[16], b, l;

    y[limbs - 1] = mpr_log_initial(x[limbs - 1]);
    for (i = 0; i < limbs - 1; i++)
        y[i] = 0;

    n = 0;
    bit_steps[n] = tol_bits;

    while (bit_steps[n] / 3 > INITIAL_PREC)
    {
        bit_steps[n+1] = bit_steps[n] / 3;
        n++;
    }

    /* Halley update for log(x): y = y + 2 - 4*t/(x + t)  where t = exp(y) */
    for ( ; n >= 0; n--)
    {
        b = bit_steps[n];
        l = (b + FLINT_BITS - 1) / FLINT_BITS;

        /* t = exp(y) */
        mpr_exp_basecase(t, y + limbs - l, l, b);

        /* u = 1 + x + t */
        u[l] = 1UL + t[l] + mpn_add_n(u, t, x + limbs - l, l);

        /* t = 4 * t, upscaled */
        mpn_lshift(t + l, t, l + 1, 2);
        mpn_zero(t, l);

        /* w = t / u */
        mpn_tdiv_q(w, t, 2*l + 1, u, l + 1);

        /* y = y + 2 - w [TODO: check for possible wraparound] */
        mpn_sub_n(y + limbs - l, y + limbs - l, w, l);
    }
}
