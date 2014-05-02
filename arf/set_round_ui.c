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

/* Top-aligns the nonzero limb v and rounds it to prec bit. */
#define ARF_NORMALISE_ROUND_LIMB(inexact, exp, v,               \
                                    sgnbit, prec, rnd)          \
    do {                                                        \
        count_leading_zeros(exp, v);                            \
        v <<= exp;                                              \
        exp = FLINT_BITS - exp;                                 \
        if (prec >= exp)                                        \
        {                                                       \
            inexact = 0;                                        \
        }                                                       \
        else                                                    \
        {                                                       \
            mp_limb_t __t = v;                                  \
            v = MASK_LIMB(v, FLINT_BITS - prec);                \
            inexact = (__t != v);                               \
            if (inexact && arf_rounds_up(rnd, sgnbit))          \
            {                                                   \
                v += (LIMB_ONE << (FLINT_BITS - prec));         \
                if (v == 0)                                     \
                {                                               \
                    v = LIMB_TOP;                               \
                    exp++;                                      \
                }                                               \
            }                                                   \
        }                                                       \
    }                                                           \
    while (0)

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
        int exp, inexact;

        ARF_NORMALISE_ROUND_LIMB(inexact, exp, v, sgnbit, prec, rnd);

        _fmpz_demote(ARF_EXPREF(x));
        ARF_EXP(x) = exp;
        ARF_DEMOTE(x);
        ARF_XSIZE(x) = ARF_MAKE_XSIZE(1, sgnbit);
        ARF_NOPTR_D(x)[0] = v;
        return inexact;
    }
}

