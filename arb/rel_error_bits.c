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

#include "arb.h"

long
arb_rel_error_bits(const arb_t x)
{
    fmpz_t midmag, radmag;
    arf_t t;
    long result;

    if (mag_is_zero(arb_radref(x)))
        return -ARF_PREC_EXACT;
    if (arf_is_special(arb_midref(x)) || mag_is_inf(arb_radref(x)))
        return ARF_PREC_EXACT;

    fmpz_init(midmag);
    fmpz_init(radmag);

    arf_abs_bound_lt_2exp_fmpz(midmag, arb_midref(x));

    arf_init_set_mag_shallow(t, arb_radref(x)); /* no need to free */
    arf_abs_bound_lt_2exp_fmpz(radmag, t);
    fmpz_add_ui(radmag, radmag, 1);

    result = _fmpz_sub_small(radmag, midmag);

    fmpz_clear(midmag);
    fmpz_clear(radmag);

    return result;
}

