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

/* todo: improve */
void
arb_add_error_2exp(arb_t x, long c)
{
    long e;

    if (COEFF_IS_MPZ(*arb_expref(x)))
    {
        printf("arb_add_error_2exp: large exponent");
        abort();
    }

    e = fmpz_get_si(arb_expref(x));

    if (c < e)
    {
        fmpz_add_ui(arb_radref(x), arb_radref(x), 1);
    }
    else
    {
        fmpz_t t;
        fmpz_init_set_ui(t, 1);
        fmpz_mul_2exp(t, t, c - e);
        fmpz_add(arb_radref(x), arb_radref(x), t);
        fmpz_clear(t);
    }
}
