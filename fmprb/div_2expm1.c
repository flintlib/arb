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



/* TODO: fix/optimize for n > prec */
void
fmprb_div_2expm1_ui(fmprb_t y, const fmprb_t x, ulong n, long prec)
{
    if (n < FLINT_BITS)
    {
        fmprb_div_ui(y, x, (1UL << n) - 1, prec);
    }
    else if (n < 1024 + prec / 32 || n > LONG_MAX / 2)
    {
        fmprb_t t;
        fmprb_init(t);
        fmprb_one(t);
        fmpz_set_ui(fmpr_expref(fmprb_midref(t)), n);
        fmprb_sub_ui(t, t, 1, prec);
        fmprb_div(y, x, t, prec);
        fmprb_clear(t);
    }
    else
    {
        fmprb_t s, t;
        long i, b;

        fmprb_init(s);
        fmprb_init(t);

        /* x / (2^n - 1) = sum_{k>=1} x * 2^(-k*n)*/

        fmprb_mul_2exp_si(s, x, -n);
        fmprb_set(t, s);
        b = 1;

        for (i = 2; i <= prec / n + 1; i++)
        {
            fmprb_mul_2exp_si(t, t, -n);
            fmprb_add(s, s, t, prec);
            b = i;
        }

        /* error bound: sum_{k>b} x * 2^(-k*n) <= x * 2^(-b*n - (n-1)) */
        fmprb_mul_2exp_si(t, x, -b*n - (n-1));
        fmprb_abs(t, t);
        fmprb_add_error(s, t);

        fmprb_set(y, s);

        fmprb_clear(s);
        fmprb_clear(t);
    }
}
