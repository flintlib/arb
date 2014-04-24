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
arf_set_trunc_mpn(arf_t y, mp_srcptr x, mp_size_t xn,
    int sgnbit, const fmpz_t topexp, mp_bitcnt_t prec)
{
    unsigned int leading;
    mp_bitcnt_t exp, bc, val, val_bits;
    mp_size_t yn, val_limbs;
    mp_ptr yptr;
    int res;

    /* Compute the total bit length of x. */
    count_leading_zeros(leading, x[xn - 1]);
    exp = xn * FLINT_BITS - leading;

    /* Set exponent. */
    if (topexp == NULL)
        fmpz_set_ui(ARF_EXPREF(y), exp);
    else
        fmpz_sub_si_inline(ARF_EXPREF(y), topexp, leading);

    /* val = mpn_scan1(x, 0) */
    val_limbs = 0;
    while (x[val_limbs] == 0)
        val_limbs++;
    count_trailing_zeros(val_bits, x[val_limbs]);
    val = val_limbs * FLINT_BITS + val_bits;

    if (exp - val <= prec)  /* Certainly exact. */
    {
        res = 0;
    }
    else    /* Certainly inexact. */
    {
        /* val = mpn_scan1(x, exp - prec) */
        mp_limb_t t;
        mp_bitcnt_t s;

        s = exp - prec;
        val_limbs = s / FLINT_BITS;
        t = x[val_limbs] & (LIMB_ONES << (s % FLINT_BITS));
        while (t == 0)
            t = x[++val_limbs];

        count_trailing_zeros(val_bits, t);
        val = val_limbs * FLINT_BITS + val_bits;

        res = 1;
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

    /* Mask off bits from the last limb. */
    if (res != 0 && prec < yn * FLINT_BITS)
    {
        yptr[0] &= LIMB_ONES << (yn * FLINT_BITS - prec);
    }

    return res;
}

