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

void
fmprb_randtest_exact(fmprb_t x, flint_rand_t state, long prec, long mag_bits)
{
    fmpr_randtest(fmprb_midref(x), state, prec, mag_bits);
    fmpr_zero(fmprb_radref(x));
}

void
fmprb_randtest_wide(fmprb_t x, flint_rand_t state, long prec, long mag_bits)
{
    fmpr_randtest(fmprb_midref(x), state, prec, mag_bits);
    fmpr_randtest(fmprb_radref(x), state, FMPRB_RAD_PREC, mag_bits);
    fmpr_abs(fmprb_radref(x), fmprb_radref(x));
}

void
fmprb_randtest_precise(fmprb_t x, flint_rand_t state, long prec, long mag_bits)
{
    fmpr_randtest(fmprb_midref(x), state, prec, mag_bits);

    if (fmpr_is_zero(fmprb_midref(x)) || (n_randint(state, 8) == 0))
    {
        fmpr_zero(fmprb_radref(x));
    }
    else
    {
        fmpz_randtest_not_zero(fmpr_manref(fmprb_radref(x)), state, FMPRB_RAD_PREC);
        fmpz_abs(fmpr_manref(fmprb_radref(x)), fmpr_manref(fmprb_radref(x)));
        fmpz_sub_ui(fmpr_expref(fmprb_radref(x)), fmpr_expref(fmprb_midref(x)), prec + FMPRB_RAD_PREC);
    }
}

void
fmprb_randtest(fmprb_t x, flint_rand_t state, long prec, long mag_bits)
{
    switch (n_randint(state, 8))
    {
        case 0:
            fmprb_randtest_exact(x, state, prec, mag_bits);
            break;
        case 1:
            fmprb_randtest_wide(x, state, prec, mag_bits);
            break;
        default:
            fmprb_randtest_precise(x, state, prec, mag_bits);
    }
}

void
fmprb_randtest_special(fmprb_t x, flint_rand_t state, long prec, long mag_bits)
{
    fmprb_randtest(x, state, prec, mag_bits);

    if (n_randint(state, 10) == 0)
        fmpr_pos_inf(fmprb_radref(x));

    switch (n_randint(state, 10))
    {
        case 0:
            fmpr_pos_inf(fmprb_midref(x));
            break;
        case 1:
            fmpr_neg_inf(fmprb_midref(x));
            break;
        case 2:
            fmpr_nan(fmprb_midref(x));
            fmpr_pos_inf(fmprb_radref(x));
            break;
        default:
            break;
    }
}

