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
_arf_get_integer_mpn(mp_ptr y, mp_srcptr x, mp_size_t xn, long exp)
{
    long bot_exp = exp - xn * FLINT_BITS;

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

void
arf_get_fmpz(fmpz_t z, const arf_t x, arf_rnd_t rnd)
{
    if (arf_is_special(x))
    {
        if (arf_is_zero(x))
        {
            fmpz_zero(z);
        }
        else
        {
            printf("arf_get_fmpz: cannot convert infinity or nan to integer\n");
            abort();
        }
    }
    else if (COEFF_IS_MPZ(*ARF_EXPREF(x)))
    {
        /* tiny */
        if (fmpz_sgn(ARF_EXPREF(x)) < 0)
        {
            int negative = ARF_SGNBIT(x);

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
        }
        else
        {
            printf("arf_get_fmpz: number too large to convert to integer\n");
            abort();
        }
    }
    else
    {
        long exp;
        int negative, inexact;
        mp_size_t xn, zn;
        mp_srcptr xp;
        __mpz_struct * zz;

        /* TBD: implement efficiently */
        if (rnd == ARF_RND_NEAR)
        {
            fmpr_t t;
            fmpr_init(t);
            arf_get_fmpr(t, x);
            fmpr_get_fmpz(z, t, rnd);
            fmpr_clear(t);
            return;
        }

        exp = ARF_EXP(x);
        negative = ARF_SGNBIT(x);

        /* |x| < 1 */
        if (exp <= 0)
        {
            if (rnd == ARF_RND_DOWN ||
                (rnd == ARF_RND_FLOOR && !negative) ||
                (rnd == ARF_RND_CEIL && negative))
            {
                fmpz_zero(z);
            }
            else
            {
                fmpz_set_si(z, negative ? -1 : 1);
            }
            return;
        }

        ARF_GET_MPN_READONLY(xp, xn, x);

        /* |x| < 2^31 or 2^63 (must save 1 bit for rounding up!) */
        if (exp < FLINT_BITS)
        {
            mp_limb_t v, v2;

            v = xp[xn - 1];
            v2 = v >> (FLINT_BITS - exp);
            inexact = (xn > 1) || ((v2 << (FLINT_BITS - exp)) != v);

            if (inexact && rnd != ARF_RND_DOWN)
            {
                if (negative && (rnd == ARF_RND_UP || rnd == ARF_RND_FLOOR))
                    v2++;
                if (!negative && (rnd == ARF_RND_UP || rnd == ARF_RND_CEIL))
                    v2++;
            }

            if (negative)
                fmpz_neg_ui(z, v2);
            else
                fmpz_set_ui(z, v2);

            return;
        }

        /* |x| >= 1 */
        zn = (exp + FLINT_BITS - 1) / FLINT_BITS;
        zz = _fmpz_promote(z);

        if (zz->_mp_alloc < zn)
            mpz_realloc2(zz, zn * FLINT_BITS);

        inexact = _arf_get_integer_mpn(zz->_mp_d, xp, xn, exp);

        zz->_mp_size = negative ? -zn : zn;
        _fmpz_demote_val(z);

        if (inexact && rnd != ARF_RND_DOWN)
        {
            if (negative && (rnd == ARF_RND_UP || rnd == ARF_RND_FLOOR))
                fmpz_sub_ui(z, z, 1);

            if (!negative && (rnd == ARF_RND_UP || rnd == ARF_RND_CEIL))
                fmpz_add_ui(z, z, 1);
        }
    }
}

