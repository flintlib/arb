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
arb_randtest(arb_t x, flint_rand_t state, long exp_bits)
{
    fmpz_randtest(arb_midref(x), state, arb_prec(x));

    if (n_randint(state, 4) == 0)
        fmpz_zero(arb_radref(x));
    else
    {
        fmpz_randtest(arb_radref(x), state, arb_prec(x));
        fmpz_abs(arb_radref(x), arb_radref(x));
    }

    if (n_randint(state, 4) == 0)
        fmpz_zero(arb_expref(x));
    else
        fmpz_randtest(arb_expref(x), state, exp_bits);

    fmpz_sub_ui(arb_expref(x), arb_expref(x), fmpz_bits(arb_midref(x)));
}
