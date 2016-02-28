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

int
_arf_set_mpn_fixed(arf_t z, mp_srcptr xp, mp_size_t xn,
        mp_size_t fixn, int negative, slong prec, arf_rnd_t rnd)
{
    slong exp, exp_shift;
    int inexact;

    exp = (slong)(xn - fixn) * FLINT_BITS;

    while (xn > 0 && xp[xn-1] == 0)
    {
        xn--;
        exp -= FLINT_BITS;
    }

    if (xn == 0)
    {
        arf_zero(z);
        return 0;
    }
    else
    {
        inexact = _arf_set_round_mpn(z, &exp_shift, xp, xn, negative, prec, rnd);
        fmpz_set_si(ARF_EXPREF(z), exp + exp_shift);
        return inexact;
    }
}

