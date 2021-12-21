/*
    Copyright (C) 2021 Albin Ahlbäck

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arf.h"

int
arf_sqr_rnd_down(arf_ptr res, arf_srcptr x, slong prec)
{
    mp_size_t xn, resn;
    mp_limb_t hi, lo;
    slong expfix;
    int sgnbit, ret, fix;
    mp_ptr resptr;

    xn = ARF_XSIZE(x);
    xn >>= 1;

    /* Either operand is a special value. */
    if (xn == 0)
    {
        arf_sqr_special(res, x);
        return 0;
    }

    if (xn == 1)
    {
        lo = ARF_NOPTR_D(x)[0];

        umul_ppmm(hi, lo, lo, lo);

        /* Shift so that the top bit is set (branch free version). */
        fix = !(hi >> (FLINT_BITS - 1));
        hi = (hi << fix) | ((lo >> (FLINT_BITS - 1)) & fix);
        lo = (lo << fix);

        ARF_DEMOTE(res);

        if (lo == 0)
        {
            resn = 1;

            if (prec >= FLINT_BITS)
            {
                lo = hi;
                ret = 0;
            }
            else
            {
                lo = MASK_LIMB(hi, FLINT_BITS - prec);
                ret = (lo != hi);
            }
        }
        else
        {
            resn = 2;

            if (prec <= FLINT_BITS)  /* Must be inexact. */
            {
                lo = MASK_LIMB(hi, FLINT_BITS - prec);
                resn = ret = 1;
            }
            else if (prec >= 2 * FLINT_BITS) /* Must be exact. */
            {
                ret = 0;
            }
            else  /* FLINT_BITS < prec < 2 * FLINT_BITS */
            {
                ret = MASK_LIMB(lo, 2 * FLINT_BITS - prec) != lo;
                lo = MASK_LIMB(lo, 2 * FLINT_BITS - prec);

                if (lo == 0)
                {
                    resn = 1;
                    lo = hi;
                }
            }
        }

        _fmpz_2times_fast(ARF_EXPREF(res), ARF_EXPREF(x), -fix);
        ARF_XSIZE(res) = ARF_MAKE_XSIZE(resn, 0);

        resptr = ARF_NOPTR_D(res);
        resptr[0] = lo;
        resptr[1] = hi;

        return ret;
    }
    else if (xn == 2)
    {
        mp_limb_t zz[4];
        mp_limb_t x1, x0;

        x0 = ARF_NOPTR_D(x)[0];
        x1 = ARF_NOPTR_D(x)[1];

        nn_sqr_2(zz[3], zz[2], zz[1], zz[0], x1, x0);

        /* Likely case, must be inexact */
        if (prec <= 2 * FLINT_BITS)
        {
            ARF_DEMOTE(res);

            fix = !(zz[3] >> (FLINT_BITS - 1));
            zz[3] = (zz[3] << fix) | ((zz[2] >> (FLINT_BITS - 1)) & fix);
            zz[2] = (zz[2] << fix) | ((zz[1] >> (FLINT_BITS - 1)) & fix);

            _fmpz_2times_fast(ARF_EXPREF(res), ARF_EXPREF(x), -fix);

            /* Rounding */
            if (prec != 2 * FLINT_BITS)
            {
                if (prec > FLINT_BITS)
                {
                    zz[2] &= (LIMB_ONES << (2 * FLINT_BITS - prec));
                }
                else if (prec == FLINT_BITS)
                {
                    zz[2] = 0;
                }
                else
                {
                    zz[3] &= (LIMB_ONES << (FLINT_BITS - prec));
                    zz[2] = 0;
                }
            }

            if (zz[2] == 0)
            {
                ARF_XSIZE(res) = ARF_MAKE_XSIZE(1, 0);
                ARF_NOPTR_D(res)[0] = zz[3];
            }
            else
            {
                ARF_XSIZE(res) = ARF_MAKE_XSIZE(2, 0);
                ARF_NOPTR_D(res)[0] = zz[2];
                ARF_NOPTR_D(res)[1] = zz[3];
            }

            return 1;
        }

        resn = 2 * xn;
        ret = _arf_set_round_mpn(res, &expfix, zz, resn, 0, prec, ARF_RND_DOWN);
        _fmpz_2times_fast(ARF_EXPREF(res), ARF_EXPREF(x), expfix);
        return ret;
    }
    else if (xn > MUL_MPFR_MIN_LIMBS && prec != ARF_PREC_EXACT
                && xn > 0.675 * prec / FLINT_BITS
                && xn < MUL_MPFR_MAX_LIMBS)  /* FIXME: proper cutoffs */
    {
        return arf_sqr_via_mpfr(res, x, prec, ARF_RND_DOWN);
    }
    else
    {
        mp_size_t alloc;
        mp_srcptr xptr;
        mp_ptr tmp;
        ARF_MUL_TMP_DECL

        ARF_GET_MPN_READONLY(xptr, xn, x);

        alloc = resn = 2 * xn;
        ARF_MUL_TMP_ALLOC(tmp, alloc)

        mpn_sqr(tmp, xptr, xn);

        ret = _arf_set_round_mpn(res, &expfix, tmp, resn, sgnbit, prec, ARF_RND_DOWN);
        _fmpz_2times_fast(ARF_EXPREF(res), ARF_EXPREF(x), expfix);
        ARF_MUL_TMP_FREE(tmp, alloc)

        return ret;
    }
}
