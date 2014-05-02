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
_arf_set_mpn_roundup_power_two(arf_t y, int sgnbit)
{
    /* The exponent has already been set, so just increment it. */
    fmpz_add_ui(ARF_EXPREF(y), ARF_EXPREF(y), 1);
    ARF_DEMOTE(y);
    ARF_NOPTR_D(y)[0] = LIMB_TOP;
    ARF_XSIZE(y) = ARF_MAKE_XSIZE(1, sgnbit);
}

int
_arf_set_round_mpn(arf_t y, long * exp_shift, mp_srcptr x, mp_size_t xn,
    int sgnbit, long prec, arf_rnd_t rnd)
{
    unsigned int leading;
    mp_bitcnt_t exp, bc, val, val_bits;
    mp_size_t yn, val_limbs;
    mp_ptr yptr;
    mp_limb_t t;
    int increment, inexact;

    /* Compute the total bit length of x. */
    count_leading_zeros(leading, x[xn - 1]);
    exp = xn * FLINT_BITS - leading;

    /* Set exponent. */
    *exp_shift = -(long) leading;

    /* Find first nonzero bit. */
    val_limbs = 0;
    while (x[val_limbs] == 0)
        val_limbs++;
    count_trailing_zeros(val_bits, x[val_limbs]);
    val = val_limbs * FLINT_BITS + val_bits;

    if (exp - val <= prec)
    {
        inexact = 0;
        increment = 0;
    }
    else
    {
        inexact = 1;
        increment = arf_rounds_up(rnd, sgnbit);

        /* Limb and bit of the truncation point. */
        val_limbs = (exp - prec) / FLINT_BITS;
        val_bits = (exp - prec) % FLINT_BITS;

        if (!increment)
        {
            /* Find first nonzero bit from the truncation point. */
            t = x[val_limbs] & (LIMB_ONES << val_bits);
            while (t == 0)
            {
                val_limbs++;
                t = x[val_limbs];
            }

            count_trailing_zeros(val_bits, t);
            val = val_limbs * FLINT_BITS + val_bits;
        }
        else
        {
            /* Find first zero bit from the truncation point */
            t = (~x[val_limbs]) & (LIMB_ONES << val_bits);
            while (t == 0)
            {
                val_limbs++;
                if (val_limbs < xn)
                    t = ~x[val_limbs];
                else  /* The array is all ones up to the highest limb. */
                {
                    val_bits = 0;
                    goto END_SCAN1;
                }
            }
            count_trailing_zeros(val_bits, t);
            END_SCAN1:
            val = val_limbs * FLINT_BITS + val_bits;

            /* Overflow to next power of two (unlikely). */
            if (val == exp)
            {
                exp_shift[0]++;
                ARF_DEMOTE(y);
                ARF_NOPTR_D(y)[0] = LIMB_TOP;
                ARF_XSIZE(y) = ARF_MAKE_XSIZE(1, sgnbit);
                return 1;
            }
        }
    }

    /* Now copy the result to destination. */
    x += val_limbs;
    xn -= val_limbs;
    bc = exp - val;
    yn = (bc + FLINT_BITS - 1) / FLINT_BITS;

    ARF_GET_MPN_WRITE(yptr, yn, y);
    ARF_XSIZE(y) |= sgnbit;

    if (leading == 0)
    {
        flint_mpn_copyi(yptr, x, xn);
    }
    else if (xn == yn)
    {
        mpn_lshift(yptr, x, yn, leading);
    }
    else
    {
        mpn_lshift(yptr, x + 1, yn, leading);
        yptr[0] |= (x[0] >> (FLINT_BITS - leading));
    }

    if (increment)
    {
        /* Mask off bits from the last limb. */
        yptr[0] &= LIMB_ONES << (yn * FLINT_BITS - bc);

        /* Increment (no carry propagation). */
        yptr[0] += LIMB_ONE << (yn * FLINT_BITS - bc);
    }
    else if (inexact && prec < yn * FLINT_BITS)
    {
        /* Mask off bits from the last limb. */
        yptr[0] &= LIMB_ONES << (yn * FLINT_BITS - prec);
    }

    return inexact;
}

