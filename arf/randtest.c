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

    Copyright (C) 2014 Fredrik Johansson

******************************************************************************/

#include "arf.h"

void
arf_randtest(arf_t x, flint_rand_t state, long bits, long mag_bits)
{
    fmpz_t t;
    fmpz_init(t);
    fmpz_randtest(t, state, bits);
    arf_set_fmpz(x, t);
    if (!arf_is_zero(x))
        fmpz_randtest(ARF_EXPREF(x), state, mag_bits);
    fmpz_clear(t);
}

void
arf_randtest_not_zero(arf_t x, flint_rand_t state, long bits, long mag_bits)
{
    fmpz_t t;
    fmpz_init(t);
    fmpz_randtest_not_zero(t, state, bits);
    arf_set_fmpz(x, t);
    fmpz_randtest(ARF_EXPREF(x), state, mag_bits);
    fmpz_clear(t);
}

void
arf_randtest_special(arf_t x, flint_rand_t state, long bits, long mag_bits)
{
    switch (n_randint(state, 32))
    {
        case 0:
            arf_zero(x);
            break;
        case 1:
            arf_pos_inf(x);
            break;
        case 2:
            arf_neg_inf(x);
            break;
        case 3:
            arf_nan(x);
            break;
        default:
            arf_randtest_not_zero(x, state, bits, mag_bits);
    }
}

