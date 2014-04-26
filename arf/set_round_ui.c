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
_arf_set_round_ui(arf_t x, ulong v, int sgnbit, long prec, arf_rnd_t rnd)
{
    if (v == 0)
    {
        arf_zero(x);
        return 0;
    }
    else
    {
        unsigned int bc, mask_bits;
        int inexact;
        mp_limb_t t;

        count_leading_zeros(bc, v);
        v <<= bc;
        bc = FLINT_BITS - bc;

        if (prec >= bc)
        {
            inexact = 0;
        }
        else
        {
            mask_bits = FLINT_BITS - prec;

            t = v;
            v = (v >> mask_bits) << mask_bits;
            inexact = (t != v);

            if (inexact && arf_rounds_up(rnd, sgnbit))
            {
                v += (LIMB_ONE << (FLINT_BITS - prec));

                /* Overflow to the next power of two (unlikely). */
                if (v == 0)
                {
                    v = LIMB_TOP;
                    bc++;
                }
            }
        }

        fmpz_set_ui(ARF_EXPREF(x), bc);
        ARF_DEMOTE(x);
        ARF_XSIZE(x) = ARF_MAKE_XSIZE(1, sgnbit);
        ARF_NOPTR_D(x)[0] = v;
        return inexact;
    }
}

