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

#include "fmpr.h"

void
fmpr_randtest(fmpr_t x, flint_rand_t state, long bits, long mag_bits)
{
    fmpz_randtest(fmpr_manref(x), state, bits);
    fmpz_randtest(fmpr_expref(x), state, mag_bits);
    fmpz_sub_ui(fmpr_expref(x), fmpr_expref(x), fmpz_bits(fmpr_manref(x)));
    _fmpr_normalise(fmpr_manref(x), fmpr_expref(x), bits, FMPR_RND_DOWN);
}

void
fmpr_randtest_not_zero(fmpr_t x, flint_rand_t state, long bits, long mag_bits)
{
    fmpz_randtest_not_zero(fmpr_manref(x), state, bits);
    fmpz_randtest(fmpr_expref(x), state, mag_bits);
    fmpz_sub_ui(fmpr_expref(x), fmpr_expref(x), fmpz_bits(fmpr_manref(x)));
    _fmpr_normalise(fmpr_manref(x), fmpr_expref(x), bits, FMPR_RND_DOWN);
}

void
fmpr_randtest_special(fmpr_t x, flint_rand_t state, long bits, long mag_bits)
{
    switch (n_randint(state, 32))
    {
        case 0:
            fmpr_zero(x);
            break;
        case 1:
            fmpr_pos_inf(x);
            break;
        case 2:
            fmpr_neg_inf(x);
            break;
        case 3:
            fmpr_nan(x);
            break;
        default:
            fmpr_randtest_not_zero(x, state, bits, mag_bits);
    }
}
