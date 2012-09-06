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

#include "fmprb.h"

/*
Fibonacci numbers using binary powering.

D. Takahashi, "A fast algorithm for computing large Fibonacci numbers",
Information Processing Letters 75 (2000) 243–246
*/

void fmprb_fib_fmpz(fmprb_t f, const fmpz_t n, long prec)
{
    fmprb_t t, u;
    long wp, sign, i;

    if (fmpz_sgn(n) < 0)
    {
        fmpz_t m;
        fmpz_init(m);
        fmpz_neg(m, n);
        fmprb_fib_fmpz(f, m, prec);
        if (fmpz_is_even(m))
            fmprb_neg(f, f);
        fmpz_clear(m);
        return;
    }

    if (fmpz_cmp_ui(n, 4) <= 0)
    {
        ulong x = fmpz_get_ui(n);
        fmprb_set_ui(f, x - (x > 1));
        return;
    }

    wp = FMPR_PREC_ADD(prec, 3 * fmpz_bits(n));

    fmprb_init(u);
    fmprb_init(t);
    fmprb_set_ui(f, 1UL);
    fmprb_set_ui(u, 1UL);
    sign = -1;

    for (i = fmpz_flog_ui(n, 2UL) - 1; i > 0; i--)
    {
        fmprb_mul(t, f, f, wp);
        fmprb_add(f, f, u, wp);
        fmprb_mul_2exp_si(f, f, -1);
        fmprb_mul(f, f, f, wp);
        fmprb_mul_2exp_si(f, f, 1);
        fmprb_submul_ui(f, t, 3, wp);
        fmprb_sub_si(f, f, 2 * sign, wp);
        fmprb_mul_ui(u, t, 5, wp);
        fmprb_add_si(u, u, 2 * sign, wp);
        sign = 1;

        if (fmpz_tstbit(n, i))
        {
            fmprb_set(t, f);
            fmprb_add(f, f, u, wp);
            fmprb_mul_2exp_si(f, f, -1);
            fmprb_mul_2exp_si(t, t, 1);
            fmprb_add(u, f, t, wp);
            sign = -1;
        }
    }

    if (fmpz_tstbit(n, 0))
    {
        fmprb_add(f, f, u, wp);
        fmprb_mul_2exp_si(f, f, -1);
        fmprb_mul(f, f, u, wp);
        fmprb_sub_si(f, f, sign, prec);
    }
    else
    {
        fmprb_mul(f, f, u, prec);
    }

    fmprb_clear(u);
    fmprb_clear(t);
}

void fmprb_fib_ui(fmprb_t f, ulong n, long prec)
{
    fmpz_t t;
    fmpz_init_set_ui(t, n);
    fmprb_fib_fmpz(f, t, prec);
    fmpz_clear(t);
}
