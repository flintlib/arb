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

#include "arb.h"

void
arb_sqrt_fmpz(arb_t x, const fmpz_t n)
{
    long size, prec;

    if (fmpz_sgn(n) < 0)
    {
        printf("arb_sqrt_fmpz: negative input\n");
        abort();
    }

    if (fmpz_is_zero(n) || fmpz_is_one(n))
    {
        arb_set_fmpz(x, n);
        return;
    }

    prec = arb_prec(x);
    size = fmpz_bits(n);

    /* todo: shift down */
    if (size >= 2 * prec)
    {
        fmpz_sqrt(arb_midref(x), n);
        fmpz_zero(arb_expref(x));
    }
    else
    {
        long shift;

        shift = 2 * prec - size;
        shift += (shift & 1);

        fmpz_mul_2exp(arb_midref(x), n, shift);
        fmpz_sqrt(arb_midref(x), arb_midref(x));
        fmpz_set_si(arb_expref(x), -shift / 2);
    }

    fmpz_one(arb_radref(x));
}
