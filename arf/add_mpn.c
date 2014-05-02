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

TLS_PREFIX mp_ptr __arf_add_tmp = NULL;
TLS_PREFIX long __arf_add_alloc = 0;

void _arf_add_tmp_cleanup(void)
{
    flint_free(__arf_add_tmp);
    __arf_add_tmp = NULL;
    __arf_add_alloc = 0;
}

/* Assumptions: top limbs of x and y nonzero. */
int
_arf_add_mpn(arf_t z, mp_srcptr xp, mp_size_t xn, int xsgnbit, const fmpz_t xexp,
                      mp_srcptr yp, mp_size_t yn, int ysgnbit, mp_bitcnt_t shift,
                      long prec, arf_rnd_t rnd)
{
    mp_size_t wn, zn, zn_original, alloc, xbase, wbase;
    mp_size_t shift_limbs;
    mp_bitcnt_t shift_bits;
    int inexact;
    long fix;
    mp_limb_t cy;
    mp_ptr tmp;
    ARF_ADD_TMP_DECL

    /* x +/- eps */
    if (shift > prec + FLINT_BITS + 1 &&
        prec != ARF_PREC_EXACT &&
        shift > (xn + 1) * FLINT_BITS)
    {
        zn = (prec + FLINT_BITS - 1) / FLINT_BITS;
        zn_original = zn = FLINT_MAX(zn, xn) + 2;
        shift_limbs = zn - xn;
        alloc = zn + 1;

        ARF_ADD_TMP_ALLOC(tmp, alloc)

        flint_mpn_zero(tmp, shift_limbs);
        flint_mpn_copyi(tmp + shift_limbs, xp, xn);

        if (xsgnbit == ysgnbit)
        {
            tmp[0] = 1;
        }
        else
        {
            mpn_sub_1(tmp, tmp, zn, 1);
            while (tmp[zn-1] == 0)
                zn--;
        }

        _arf_set_round_mpn(z, &fix, tmp, zn, xsgnbit, prec, rnd);
        fix += (zn - zn_original) * FLINT_BITS;
        _fmpz_add_fast(ARF_EXPREF(z), xexp, fix);

        ARF_ADD_TMP_FREE(tmp, alloc)
        return 1;
    }

    shift_limbs = shift / FLINT_BITS;
    shift_bits = shift % FLINT_BITS;

    /* TODO: shift == 0 optimization! (common case) */

    /* w = y shifted */
    wn = yn + (shift_bits != 0);

    /* Size of sum (possibly one more limb than this for carry). */
    zn_original = zn = FLINT_MAX(xn, wn + shift_limbs);
    alloc = zn + 1;

    ARF_ADD_TMP_ALLOC(tmp, alloc)

    /* Low limbs of both operands */
    xbase = zn - xn;
    wbase = zn - shift_limbs - wn;

    /* Shift y to its correct offset in the sum. */
    if (shift_bits == 0)
        flint_mpn_copyi(tmp + wbase, yp, yn);
    else
        tmp[wbase] = cy = mpn_rshift(tmp + wbase + 1, yp, yn, shift_bits);

    if (xsgnbit == ysgnbit)
    {
        if (shift_limbs >= xn)
        {
            /*   XXXXXX
                         WWWWWWW    */
            flint_mpn_zero(tmp + wn, zn - xn - wn);
            flint_mpn_copyi(tmp + zn - xn, xp, xn);
            cy = 0;
        }
        else if (shift_limbs == 0)
        {
            /*  XXXXXX        or    XXXXX
                WWWW                WWWWWWW     */
            if (wn < xn)
            {
                flint_mpn_copyi(tmp, xp, xn - wn);
                cy = mpn_add_n(tmp + wbase, xp + wbase, tmp + wbase, wn);
            }
            else
            {
                cy = mpn_add_n(tmp + xbase, xp, tmp + xbase, xn);
            }
        }
        else if (xbase >= wbase)
        {
            /* XXXXXX         or  XXXXXX
                 WWWWWWW            WWWW    */
            cy = mpn_add_n(tmp + xbase, tmp + xbase, xp, wn - xbase);
            cy = mpn_add_1(tmp + wn, xp + wn - xbase, zn - wn, cy);
        }
        else
        {
            /* XXXXXXXX           todo: merge with above?
                  WWW    */
            flint_mpn_copyi(tmp, xp, wbase);
            cy = mpn_add_n(tmp + wbase, tmp + wbase, xp + wbase, wn);
            cy = mpn_add_1(tmp + zn - shift_limbs, xp + xn - shift_limbs, shift_limbs, cy);
        }

        /* There could be a carry. */
        tmp[zn] = cy;
        zn += cy;
    }
    else
    {
        if (shift_limbs >= xn)
        {
            /*   XXXXXX
                         WWWWWWW    */
            mpn_neg(tmp, tmp, wn);
            flint_mpn_store(tmp + wn, zn - xn - wn, -LIMB_ONE);
            mpn_sub_1(tmp + zn - xn, xp, xn, 1);
        }
        else if (shift_limbs == 0)
        {
            /*  XXXXXX (1)    or    XXXXX   (2)
                WWWW                WWWWWWW     */
            if (wn <= xn) /* (1) */
            {
                if (mpn_cmp(xp + wbase, tmp + wbase, wn) >= 0)
                {
                    flint_mpn_copyi(tmp, xp, wbase);
                    mpn_sub_n(tmp + wbase, xp + wbase, tmp + wbase, wn);
                }
                else
                {
                    xsgnbit = ysgnbit;
                    mpn_sub_n(tmp + wbase, tmp + wbase, xp + wbase, wn);
                    if (wbase != 0)
                    {
                        cy = mpn_neg_n(tmp, xp, wbase);
                        mpn_sub_1(tmp + wbase, tmp + wbase, wn, cy);
                    }
                }
            }
            else /* (2) */
            {
                if (mpn_cmp(tmp + xbase, xp, xn) >= 0)
                {
                    xsgnbit = ysgnbit;
                    mpn_sub_n(tmp + xbase, tmp + xbase, xp, xn);
                }
                else
                {
                    cy = mpn_neg(tmp, tmp, xbase);
                    mpn_sub_n(tmp + xbase, xp, tmp + xbase, xn);
                    mpn_sub_1(tmp + xbase, tmp + xbase, xn, cy);
                }
            }
        }
        else if (xbase >= wbase)
        {
            /* XXXXXX         or  XXXXXX
                 WWWWWWW            WWWW    */
            if (xbase > 0)
            {
                cy = mpn_sub_n(tmp + xbase, xp, tmp + xbase, wn - xbase);
                cy = mpn_sub_1(tmp + wn, xp + xn - shift_limbs, shift_limbs, cy);
                cy = mpn_neg(tmp, tmp, xbase);
                mpn_sub_1(tmp + xbase, tmp + xbase, xn, cy);
            }
            else
            {
                cy = mpn_sub_n(tmp, xp, tmp, wn);
                cy = mpn_sub_1(tmp + wn, xp + wn, xn - wn, cy);
            }
        }
        else
        {
            /* XXXXXXXX
                  WWW    */
            flint_mpn_copyi(tmp, xp, wbase);
            cy = mpn_sub_n(tmp + wbase, xp + wbase, tmp + wbase, wn);
            mpn_sub_1(tmp + zn - shift_limbs, xp + xn - shift_limbs, shift_limbs, cy);
        }

        /* There could be cancellation. */
        while (zn > 0 && tmp[zn - 1] == 0)
            zn--;
    }

    if (zn == 0)
    {
        arf_zero(z);
        inexact = 0;
    }
    else
    {
        inexact = _arf_set_round_mpn(z, &fix, tmp, zn, xsgnbit, prec, rnd);
        fix += (zn - zn_original) * FLINT_BITS;
        _fmpz_add_fast(ARF_EXPREF(z), xexp, fix);
    }

    ARF_ADD_TMP_FREE(tmp, alloc)
    return inexact;
}

