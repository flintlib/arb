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
Computes the exact interval represented by x, in the form of an integer
interval multiplied by a power of two, i.e. x = [a, b] * 2^exp.

The outcome is undefined if the midpoint or radius of x is non-finite,
or if the difference in magnitude between the midpoint and radius
is so large that representing the endpoints exactly would cause overflows.
*/

void
fmprb_get_interval_fmpz_2exp(fmpz_t a, fmpz_t b, fmpz_t exp, const fmprb_t x)
{

    if (fmprb_is_exact(x))
    {
        fmpr_get_fmpz_2exp(a, exp, fmprb_midref(x));
        fmpz_set(b, a);
    }
    else
    {
        fmpr_t t;
        fmpz_t exp2;
        long s;

        fmpr_init(t);
        fmpz_init(exp2);

        fmpr_sub(t, fmprb_midref(x), fmprb_radref(x), FMPR_PREC_EXACT, FMPR_RND_DOWN);
        fmpr_get_fmpz_2exp(a, exp, t);

        fmpr_add(t, fmprb_midref(x), fmprb_radref(x), FMPR_PREC_EXACT, FMPR_RND_DOWN);
        fmpr_get_fmpz_2exp(b, exp2, t);

        s = _fmpz_sub_small(exp, exp2);

        if (s <= 0)
        {
            fmpz_mul_2exp(b, b, -s);
        }
        else
        {
            fmpz_mul_2exp(a, a, s);
            fmpz_set(exp, exp2);
        }

        fmpr_clear(t);
        fmpz_clear(exp2);
    }
}
