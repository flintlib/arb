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

/* assumes xp[0] is nonzero and xp[xn-1] has the top bit set */
/* assumes that the correct number of limbs for y has been allocated
   (no zero-extension is done) */
/* truncates, returns inexact */
int
_arf_get_integer_mpn(mp_ptr y, mp_srcptr x, mp_size_t xn, slong exp)
{
    slong bot_exp = exp - xn * FLINT_BITS;

    if (bot_exp >= 0)
    {
        mp_size_t bot_limbs;
        mp_bitcnt_t bot_bits;

        bot_limbs = bot_exp / FLINT_BITS;
        bot_bits = bot_exp % FLINT_BITS;

        flint_mpn_zero(y, bot_limbs);

        if (bot_bits == 0)
            flint_mpn_copyi(y + bot_limbs, x, xn);
        else
            y[bot_limbs + xn] = mpn_lshift(y + bot_limbs, x, xn, bot_bits);

        /* exact */
        return 0;
    }
    else if (exp <= 0)
    {
        /* inexact */
        return 1;
    }
    else
    {
        mp_size_t top_limbs;
        mp_bitcnt_t top_bits;
        mp_limb_t cy;

        top_limbs = exp / FLINT_BITS;
        top_bits = exp % FLINT_BITS;

        if (top_bits == 0)
        {
            flint_mpn_copyi(y, x + xn - top_limbs, top_limbs);
            /* inexact */
            return 1;
        }
        else
        {
            /* can be inexact */
            cy = mpn_rshift(y, x + xn - top_limbs - 1,
                top_limbs + 1, FLINT_BITS - top_bits);

            return (cy != 0) || (top_limbs + 1 != xn);
        }
    }
}

int
arf_get_fmpz(fmpz_t z, const arf_t x, arf_rnd_t rnd)
{
    slong exp;
    int negative, inexact, value, roundup;
    mp_size_t xn, zn;
    mp_srcptr xp;
    __mpz_struct * zz;
    mp_ptr zp;
    mp_limb_t v, v2, v3;

    if (arf_is_special(x))
    {
        if (arf_is_zero(x))
        {
            fmpz_zero(z);
            return 0;
        }
        else
        {
            flint_printf("arf_get_fmpz: cannot convert infinity or nan to integer\n");
            abort();
        }
    }

    exp = ARF_EXP(x);
    negative = ARF_SGNBIT(x);

    if (COEFF_IS_MPZ(exp))
    {
        /* tiny */
        if (fmpz_sgn(ARF_EXPREF(x)) < 0)
        {
            if (rnd == ARF_RND_NEAR
                || rnd == ARF_RND_DOWN
                || (rnd == ARF_RND_FLOOR && !negative)
                || (rnd == ARF_RND_CEIL && negative))
            {
                fmpz_zero(z);
            }
            else
            {
                fmpz_set_si(z, negative ? -1 : 1);
            }

            return 1;
        }
        else
        {
            flint_printf("arf_get_fmpz: number too large to convert to integer\n");
            abort();
        }
    }

    /* |x| < 1 */
    if (exp <= 0)
    {
        if (rnd == ARF_RND_NEAR)
        {
            if (exp == 0)
            {
                /* check for the special case +/- 1/2 */
                ARF_GET_MPN_READONLY(xp, xn, x);

                if (xp[xn - 1] < LIMB_TOP || (xn == 1 && xp[xn - 1] == LIMB_TOP))
                    value = 0;
                else
                    value = negative ? -1 : 1;
            }
            else
            {
                value = 0;
            }
        }
        else if (rnd == ARF_RND_DOWN ||
            (rnd == ARF_RND_FLOOR && !negative) ||
            (rnd == ARF_RND_CEIL && negative))
        {
            value = 0;
        }
        else
        {
            value = negative ? -1 : 1;
        }

        _fmpz_demote(z);
        *z = value;
        return 1;
    }

    ARF_GET_MPN_READONLY(xp, xn, x);

    /* Fast case: |x| < 2^31 or 2^63 (must save 1 bit for rounding up!) */
    if (exp < FLINT_BITS)
    {
        v = xp[xn - 1];
        v2 = v >> (FLINT_BITS - exp); /* integral part */
        v3 = v << exp;                /* fractional part (truncated, at least 1 bit) */
        inexact = (xn > 1) || (v3 != 0);

        if (inexact && rnd != ARF_RND_DOWN)
        {
            if (rnd == ARF_RND_NEAR)
            {
                /* round up of fractional part is > 1/2,
                   or if equal to 1/2 and the integral part is odd */
                v2 += (v3 > LIMB_TOP) || (v3 == LIMB_TOP && (xn > 1 || (v2 & 1)));
            }
            else
            {
                v2 += (rnd == ARF_RND_UP) || (negative ^ (rnd == ARF_RND_CEIL));
            }
        }

        if (negative)
            fmpz_neg_ui(z, v2);
        else
            fmpz_set_ui(z, v2);
        return inexact;
    }

    /* |x| >= 1 */

    /* Allocate space for result + 1 extra bit. We need one extra bit
       temporarily to check rounding to nearest. We also need one extra bit
       to round up. */
    zn = (exp + (rnd != ARF_RND_DOWN) + FLINT_BITS - 1) / FLINT_BITS;

    zz = _fmpz_promote(z);
    if (zz->_mp_alloc < zn)
        mpz_realloc2(zz, zn * FLINT_BITS);

    zp = zz->_mp_d;

    if (rnd == ARF_RND_DOWN)
    {
        /* zn is the exact size */
        inexact = _arf_get_integer_mpn(zp, xp, xn, exp);
    }
    else
    {
        zp[zn - 1] = 0;
        inexact = _arf_get_integer_mpn(zp, xp, xn, exp + (rnd == ARF_RND_NEAR));

        if (rnd == ARF_RND_NEAR)
        {
            v = zp[0];
            /* round up if fractional part is >= 1/2 and (there are
               more discarded bits, or the truncated value would be odd) */
            roundup = (v & 1) & (inexact | (v >> 1));
            inexact |= (v & 1);
            mpn_rshift(zp, zp, zn, 1);
            mpn_add_1(zp, zp, zn, roundup);
        }
        else if (inexact && ((rnd == ARF_RND_UP)
            || (negative ^ (rnd == ARF_RND_CEIL))))
        {
            mpn_add_1(zp, zp, zn, 1);
        }

        zn -= (zp[zn - 1] == 0);
    }

    zz->_mp_size = negative ? -zn : zn;
    _fmpz_demote_val(z);
    return inexact;
}

int
arf_get_fmpz_fixed_fmpz(fmpz_t y, const arf_t x, const fmpz_t e)
{
    if (arf_is_special(x))
    {
        return arf_get_fmpz(y, x, ARF_RND_DOWN);
    }
    else
    {
        int inexact;
        fmpz_t exp;
        arf_t tmp;

        fmpz_init(exp);
        fmpz_sub(exp, ARF_EXPREF(x), e);
        arf_init_set_shallow(tmp, x);
        ARF_EXP(tmp) = *exp;
        inexact = arf_get_fmpz(y, tmp, ARF_RND_DOWN);
        fmpz_clear(exp);
        return inexact;
    }
}

int
arf_get_fmpz_fixed_si(fmpz_t y, const arf_t x, slong e)
{
    if (arf_is_special(x))
    {
        return arf_get_fmpz(y, x, ARF_RND_DOWN);
    }
    else
    {
        int inexact;
        fmpz_t exp;
        arf_t tmp;

        fmpz_init(exp);
        fmpz_sub_si(exp, ARF_EXPREF(x), e);
        arf_init_set_shallow(tmp, x);
        ARF_EXP(tmp) = *exp;
        inexact = arf_get_fmpz(y, tmp, ARF_RND_DOWN);
        fmpz_clear(exp);
        return inexact;
    }
}

