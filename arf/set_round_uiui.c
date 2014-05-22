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
_arf_set_round_uiui(arf_t z, long * fix, mp_limb_t hi, mp_limb_t lo, int sgnbit, long prec, arf_rnd_t rnd)
{
    int leading, trailing, bc, inexact, zn, up, exp;

    if (hi == 0)
    {
        ARF_NORMALISE_ROUND_LIMB(inexact, exp, lo, sgnbit, prec, rnd);
        leading = 2 * FLINT_BITS - exp;
        zn = 1;
    }
    else if (lo == 0)
    {
        ARF_NORMALISE_ROUND_LIMB(inexact, exp, hi, sgnbit, prec, rnd);
        leading = FLINT_BITS - exp;
        lo = hi;
        zn = 1;
    }
    else
    {
        count_leading_zeros(leading, hi);
        count_trailing_zeros(trailing, lo);

        bc = 2 * FLINT_BITS - leading - trailing;

        if (bc <= prec)
        {
            inexact = 0;
            zn = 2;

            if (leading != 0)
            {
                if (bc <= FLINT_BITS)
                {
                    lo = (hi << leading) | (lo >> (FLINT_BITS - leading));
                    zn = 1;
                }
                else
                {
                    hi = (hi << leading) | (lo >> (FLINT_BITS - leading));
                    lo = (lo << leading);
                }
            }
        }
        else
        {
            inexact = 1;
            up = arf_rounds_up(rnd, sgnbit);

            if (prec <= FLINT_BITS)
            {
                zn = 1;

                if (leading == 0)
                    lo = hi;
                else
                    lo = (hi << leading) | (lo >> (FLINT_BITS - leading));

                lo = lo & (LIMB_ONES << (FLINT_BITS - prec));

                if (up)
                {
                    mp_limb_t t, ovf;
                    t = lo + (LIMB_ONE << (FLINT_BITS - prec));
                    ovf = (t == 0);
                    leading -= ovf;
                    lo = (t >> ovf) | (ovf << (FLINT_BITS - 1));
                }
            }
            else
            {
                zn = 2;

                if (leading != 0)
                {
                    hi = (hi << leading) | (lo >> (FLINT_BITS - leading));
                    lo = (lo << leading);
                }

                lo = lo & (LIMB_ONES << (2 * FLINT_BITS - prec));

                if (up)
                {
                    add_ssaaaa(hi, lo, hi, lo, 0, (LIMB_ONE << (2 * FLINT_BITS - prec)));
                }

                if (lo == 0)
                {
                    if (hi == 0)
                    {
                        leading -= 1;
                        lo = LIMB_TOP;
                    }
                    else
                    {
                        lo = hi;
                    }

                    zn = 1;
                }
            }
        }
    }

    *fix = -leading;
    ARF_DEMOTE(z);
    ARF_XSIZE(z) = ARF_MAKE_XSIZE(zn, sgnbit);
    ARF_NOPTR_D(z)[0] = lo;
    ARF_NOPTR_D(z)[1] = hi;
    return inexact;
}

