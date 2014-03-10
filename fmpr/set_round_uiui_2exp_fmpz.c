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

    Copyright (C) 2013 Fredrik Johansson

******************************************************************************/

#include "fmpr.h"

long
fmpr_set_round_uiui_2exp_fmpz(fmpr_t z,
    mp_limb_t hi, mp_limb_t lo, const fmpz_t exp, int negative,
    long prec, fmpr_rnd_t rnd)
{
    if (hi == 0)
    {
        return fmpr_set_round_ui_2exp_fmpz(z, lo, exp, negative, prec, rnd);
    }
    else
    {
        /* 2 limbs */
        long lead, trail, bc, shift, ret;

        if ((lo & 1) == 0)
        {
            if (lo == 0)
            {
                shift = FLINT_BITS;
                ret = fmpr_set_round_ui_2exp_fmpz(z, hi, exp, negative, prec, rnd);
            }
            else
            {
                count_trailing_zeros(shift, lo);
                lo = (hi << (FLINT_BITS - shift)) | lo >> shift;
                hi = hi >> shift;
                ret = fmpr_set_round_uiui_2exp_fmpz(z, hi, lo, exp, negative, prec, rnd);
            }

            fmpz_add_si_inline(fmpr_expref(z), fmpr_expref(z), shift);
            return ret;
        }

        count_leading_zeros(lead, hi);
        bc = 2 * FLINT_BITS - lead;

        shift = 0;
        ret = FMPR_RESULT_EXACT;

        if (bc > prec)
        {
            shift = bc - prec;

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
                ret -= (trail + FLINT_BITS == prec);

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

        fmpz_add_si_inline(fmpr_expref(z), exp, shift);
        return ret;
    }
}

