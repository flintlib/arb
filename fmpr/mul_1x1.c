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

    Copyright (C) 2013, 2014 Fredrik Johansson

******************************************************************************/

#include "fmpr.h"

long
_fmpr_mul_1x1(fmpr_t z, mp_limb_t u, const fmpz_t xexp, mp_limb_t v,
    const fmpz_t yexp, int negative, long prec, fmpr_rnd_t rnd)
{
    long lead, trail, bc, shift, ret;
    mp_limb_t hi, lo;

    umul_ppmm(hi, lo, u, v);
    shift = 0;

    if (hi == 0)
    {
        /* 1 limb */
        count_leading_zeros(lead, lo);
        bc = FLINT_BITS - lead;

        ret = FMPR_RESULT_EXACT;
        if (bc > prec)
        {
            shift += bc - prec;
            lo = (lo >> shift) + rounds_up(rnd, negative);
            count_trailing_zeros(trail, lo);
            lo >>= trail;
            shift += trail;
            ret = trail;

            /* special case: if the mantissa overflowed to the next power of two,
               the error bound must be multiplied by two */
            ret -= (trail == prec);
        }

        if (!negative)
            fmpz_set_ui(fmpr_manref(z), lo);
        else
            fmpz_neg_ui(fmpr_manref(z), lo);
    }
    else
    {
        /* 2 limbs */
        count_leading_zeros(lead, hi);
        bc = 2 * FLINT_BITS - lead;

        ret = FMPR_RESULT_EXACT;

        if (bc > prec)
        {
            shift += bc - prec;

            /* round */
            if (shift < FLINT_BITS)
            {
                lo = (lo >> shift) | (hi << (FLINT_BITS - shift));
                hi >>= shift;
            }
            else
            {
                lo = hi >> (shift - FLINT_BITS);
                hi = 0;
            }

            if (rounds_up(rnd, negative))
                add_ssaaaa(hi, lo, hi, lo, 0, 1);

            /* remove trailing zeros */
            if (lo == 0)
            {
                count_trailing_zeros(trail, hi);
                hi >>= trail;
                shift += FLINT_BITS + trail;
                ret = FLINT_BITS + trail;

                /* special case: if the mantissa overflowed to the next power of two,
                   the error bound must be multiplied by two */
                ret -= (FLINT_BITS + trail == prec);

                if (!negative)
                    fmpz_set_ui(fmpr_manref(z), hi);
                else
                    fmpz_neg_ui(fmpr_manref(z), hi);
            }
            else
            {
                count_trailing_zeros(trail, lo);

                if (trail != 0)
                {
                    lo = (lo >> trail) | (hi << (FLINT_BITS - trail));
                    hi >>= trail;
                    shift += trail;
                }
                ret = trail;

                /* special case: if the mantissa overflowed to the next power of two,
                   the error bound must be multiplied by two */
                ret -= (trail == prec);

                if (!negative)
                    fmpz_set_uiui(fmpr_manref(z), hi, lo);
                else
                    fmpz_neg_uiui(fmpr_manref(z), hi, lo);
            }
        }
        else
        {
            if (!negative)
                fmpz_set_uiui(fmpr_manref(z), hi, lo);
            else
                fmpz_neg_uiui(fmpr_manref(z), hi, lo);
        }
    }

    fmpz_add2_fmpz_si_inline(fmpr_expref(z), xexp, yexp, shift);
    return ret;
}

