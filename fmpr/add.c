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

    Copyright (C) 2012, 2013 Fredrik Johansson

******************************************************************************/

#include "fmpr.h"

#define ADD_STACK_ALLOC 40
#define ADD_TLS_ALLOC 1000

TLS_PREFIX mp_ptr __add_tmp = NULL;
TLS_PREFIX long __add_alloc = 0;

void _add_tmp_cleanup(void)
{
    flint_free(__add_tmp);
    __add_alloc = 0;
}

#define ADD_TMP_ALLOC \
    if (alloc <= ADD_STACK_ALLOC) \
    { \
        tmp = tmp_stack; \
    } \
    else if (alloc <= ADD_TLS_ALLOC) \
    { \
        if (__add_alloc < alloc) \
        { \
            if (__add_alloc == 0) \
            { \
                flint_register_cleanup_function(_add_tmp_cleanup); \
            } \
            __add_tmp = flint_realloc(__add_tmp, sizeof(mp_limb_t) * alloc); \
            __add_alloc = alloc; \
        } \
        tmp = __add_tmp; \
    } \
    else \
    { \
        tmp = flint_malloc(sizeof(mp_limb_t) * alloc); \
    }

#define ADD_TMP_FREE \
    if (alloc > ADD_TLS_ALLOC) \
        flint_free(tmp);

#define READ_SIGN_PTR(zsign, ztmp, zptr, zn, zv) \
    if (!COEFF_IS_MPZ(zv)) \
    { \
        zsign = zv < 0; \
        ztmp = FLINT_ABS(zv); \
        zptr = &ztmp; \
        zn = 1; \
    } \
    else \
    { \
        __mpz_struct * zz = COEFF_TO_PTR(zv); \
        zptr = zz->_mp_d; \
        zn = zz->_mp_size; \
        zsign = zn < 0; \
        zn = FLINT_ABS(zn); \
    }

static __inline__ long
_fmpr_add_small(fmpr_t z, mp_limb_t x, int xsign, const fmpz_t xexp, mp_limb_t y, int ysign, const fmpz_t yexp, long shift, long prec, long rnd)
{
    mp_limb_t hi, lo, t, u;
    int sign = ysign;

    t = x;
    u = y;

    lo = u << shift;
    hi = (shift == 0) ? 0 : (u >> (FLINT_BITS - shift));

    if (xsign == ysign)
    {
        add_ssaaaa(hi, lo, hi, lo, 0, t);
    }
    else
    {
        if (hi == 0)
        {
            if (lo >= t)
            {
                lo = lo - t;
            }
            else
            {
                lo = t - lo;
                sign = !sign;
            }
        }
        else
        {
            sub_ddmmss(hi, lo, hi, lo, 0, t);
        }
    }

    if (hi == 0)
        return fmpr_set_round_ui_2exp_fmpz(z, lo, xexp, sign, prec, rnd);
    else
        return fmpr_set_round_uiui_2exp_fmpz(z, hi, lo, xexp, sign, prec, rnd);
}

/* computes x + y * 2^shift (optionally negated) */
long
_fmpr_add_large(fmpr_t z,
        mp_srcptr xman, mp_size_t xn, int xsign, const fmpz_t xexp,
        mp_srcptr yman, mp_size_t yn, int ysign, const fmpz_t yexp,
        long shift, long prec, fmpr_rnd_t rnd)
{
    long tn, zn, alloc, ret, shift_bits, shift_limbs;
    int negative;
    mp_limb_t tmp_stack[ADD_STACK_ALLOC];
    mp_limb_t cy;
    mp_ptr tmp, tmp2;

    shift_limbs = shift / FLINT_BITS;
    shift_bits = shift % FLINT_BITS;

    /* x does not overlap with y or the result -- outcome is
       equivalent to adding/subtracting a small number to/from y
       and rounding */
    if (shift > xn * FLINT_BITS &&
        prec != FMPR_PREC_EXACT &&
        xn * FLINT_BITS + prec - (FLINT_BITS * (yn - 1)) < shift)
    {
        zn = (prec + FLINT_BITS - 1) / FLINT_BITS;
        zn = FLINT_MAX(zn, yn) + 2;
        shift_limbs = zn - yn;
        alloc = zn;

        ADD_TMP_ALLOC

        flint_mpn_zero(tmp, shift_limbs);
        flint_mpn_copyi(tmp + shift_limbs, yman, yn);

        if (xsign == ysign)
        {
            tmp[0] = 1;
        }
        else
        {
            mpn_sub_1(tmp, tmp, zn, 1);
            while (tmp[zn-1] == 0)
                zn--;
        }

        ret = _fmpr_set_round_mpn(&shift, fmpr_manref(z), tmp, zn, ysign, prec, rnd);
        shift -= shift_limbs * FLINT_BITS;
        fmpz_add_si_inline(fmpr_expref(z), yexp, shift);

        ADD_TMP_FREE

        return ret;
    }

    if (xsign == ysign)
    {
        negative = xsign;

        /* certainly no overlap between x and y * 2^shift

                        XXXXX
               YYYYY         
                    |-------|
                   shift_limbs

          todo: handle huge case efficiently!
        */
        if (shift_limbs >= xn)
        {
            alloc = shift_limbs + yn + 1;

            ADD_TMP_ALLOC

            flint_mpn_copyi(tmp, xman, xn);
            flint_mpn_zero(tmp + xn, shift_limbs - xn);

            if (shift_bits == 0)
            {
                flint_mpn_copyi(tmp + shift_limbs, yman, yn);
                zn = shift_limbs + yn;
            }
            else
            {
                cy = mpn_lshift(tmp + shift_limbs, yman, yn, shift_bits);
                tmp[shift_limbs + yn] = cy;
                zn = shift_limbs + yn + (cy != 0);
            }
        }
        else /* may have overlap between x and y * 2^shift */
        {
            long alloc1, alloc2;

            alloc2 = yn + 1;   /* shifted value of y */
            alloc1 = FLINT_MAX(alloc2 + shift_limbs, xn) + 1; /* space for sum */
            alloc = alloc1 + alloc2;

            ADD_TMP_ALLOC

            tmp2 = tmp + alloc1;

            flint_mpn_copyi(tmp, xman, shift_limbs);

            /* {t, tn} = shift of y (ignoring shift_limbs zero limbs) */
            if (shift_bits == 0)
            {
                flint_mpn_copyi(tmp2, yman, yn);
                tn = yn;
            }
            else
            {
                cy = mpn_lshift(tmp2, yman, yn, shift_bits);
                tmp2[yn] = cy;
                tn = yn + (cy != 0);
            }

            if (tn + shift_limbs == xn)
            {
                /*    XXXXXXXX
                      TTTTT        */
                cy = mpn_add_n(tmp + shift_limbs, xman + shift_limbs, tmp2, tn);
                tmp[xn] = cy;
                zn = xn + (cy != 0);
            }
            else if (tn + shift_limbs < xn)
            {
                /*      XXXXXXXX
                          TTTTT    */
                cy = mpn_add_n(tmp + shift_limbs, xman + shift_limbs, tmp2, tn);
                cy = mpn_add_1(tmp + shift_limbs + tn, xman + shift_limbs + tn, xn - tn - shift_limbs, cy);
                tmp[xn] = cy;
                zn = xn + (cy != 0);
            }
            else
            {
                /*    XXXXXXXX
                   TTTTTTT         */
                cy = mpn_add_n(tmp + shift_limbs, xman + shift_limbs, tmp2, xn - shift_limbs);
                cy = mpn_add_1(tmp + xn, tmp2 + xn - shift_limbs, tn + shift_limbs - xn, cy);
                tmp[tn + shift_limbs] = cy;
                zn = tn + shift_limbs + (cy != 0);
            }
        }
    }
    else
    {
        /* certainly no overlap between x and y * 2^shift

                        XXXXX
               YYYYY         
                    |-------|
                   shift_limbs

          todo: handle huge case efficiently!
        */
        if (shift_limbs >= xn)
        {
            alloc = shift_limbs + yn + 1;

            ADD_TMP_ALLOC

            mpn_neg_n(tmp, xman, xn);

            /* since x is nonzero, there is certainly carry (borrow) */
            flint_mpn_store(tmp + xn, shift_limbs - xn, (mp_limb_t) -1);

            if (shift_bits == 0)
            {
                flint_mpn_copyi(tmp + shift_limbs, yman, yn);
                zn = shift_limbs + yn;
            }
            else
            {
                cy = mpn_lshift(tmp + shift_limbs, yman, yn, shift_bits);
                tmp[shift_limbs + yn] = cy;
                zn = shift_limbs + yn + (cy != 0);
            }

            mpn_sub_1(tmp + shift_limbs, tmp + shift_limbs, zn - shift_limbs, 1);

            /* y has larger absolute value, and determines the sign */
            negative = ysign;
        }
        else /* may have overlap between x and y * 2^shift */
        {
            long alloc1, alloc2;

            alloc2 = yn + 1;   /* shifted value of y */
            alloc1 = FLINT_MAX(alloc2 + shift_limbs, xn); /* space for difference */
            alloc = alloc1 + alloc2;

            ADD_TMP_ALLOC

            tmp2 = tmp + alloc1;

            /* {t, tn} = shift of y (ignoring shift_limbs zero limbs) */
            if (shift_bits == 0)
            {
                flint_mpn_copyi(tmp2, yman, yn);
                tn = yn;
            }
            else
            {
                cy = mpn_lshift(tmp2, yman, yn, shift_bits);
                tmp2[yn] = cy;
                tn = yn + (cy != 0);
            }

            if (tn + shift_limbs == xn)
            {
                int comp;
                /*    XXXXXXXX
                      TTTTT        */

                comp = mpn_cmp(xman + shift_limbs, tmp2, tn);

                if (comp == 0)
                {
                    long i;

                    for (i = 0; i < shift_limbs; i++)
                    {
                        if (xman[i] != 0)
                        {
                            comp = 1;
                            break;
                        }
                    }

                    /* the result is exactly zero */
                    if (comp == 0)
                    {
                        fmpr_zero(z);
                        ret = FMPR_RESULT_EXACT;
                        goto cleanup;
                    }
                }

                /* x is larger */
                if (comp > 0)
                {
                    flint_mpn_copyi(tmp, xman, shift_limbs);
                    mpn_sub_n(tmp + shift_limbs, xman + shift_limbs, tmp2, tn);
                    negative = xsign;
                }
                else
                {
                    mpn_sub_n(tmp + shift_limbs, tmp2, xman + shift_limbs, tn);
                    if (shift_limbs != 0)
                    {
                        cy = mpn_neg_n(tmp, xman, shift_limbs);
                        mpn_sub_1(tmp + shift_limbs, tmp + shift_limbs, tn, cy);
                    }
                    negative = ysign;
                }

                zn = xn;
            }
            else if (tn + shift_limbs < xn)
            {
                /*      XXXXXXXX
                          TTTTT    */

                flint_mpn_copyi(tmp, xman, shift_limbs);
                cy = mpn_sub_n(tmp + shift_limbs, xman + shift_limbs, tmp2, tn);
                mpn_sub_1(tmp + shift_limbs + tn, xman + shift_limbs + tn, xn - tn - shift_limbs, cy);
                zn = xn;
                negative = xsign;
            }
            else
            {
                /*    XXXXXXXX
                   TTTTTTT         */
                cy = mpn_sub_n(tmp + shift_limbs, tmp2, xman + shift_limbs, xn - shift_limbs);
                mpn_sub_1(tmp + xn, tmp2 + xn - shift_limbs, tn + shift_limbs - xn, cy);

                if (shift_limbs != 0)
                {
                    cy = mpn_neg_n(tmp, xman, shift_limbs);
                    mpn_sub_1(tmp + shift_limbs, tmp + shift_limbs, tn, cy);
                }

                zn = shift_limbs + tn;
                negative = ysign;
            }
        }

        /* there might have been cancellation */
        while (tmp[zn-1] == 0)
            zn--;
    }

    ret = _fmpr_set_round_mpn(&shift, fmpr_manref(z), tmp, zn, negative, prec, rnd);
    fmpz_add_si_inline(fmpr_expref(z), xexp, shift);

cleanup:
    ADD_TMP_FREE

    return ret;
}

static long
_fmpr_add_special(fmpr_t z, const fmpr_t x, const fmpr_t y, long prec, fmpr_rnd_t rnd)
{
    if (fmpr_is_zero(x))
    {
        if (fmpr_is_zero(y))
        {
            fmpr_zero(z);
            return FMPR_RESULT_EXACT;
        }
        else
            return fmpr_set_round(z, y, prec, rnd);
    }
    else if (fmpr_is_zero(y))
    {
        return fmpr_set_round(z, x, prec, rnd);
    }
    else if (fmpr_is_nan(x) || fmpr_is_nan(y)
        || (fmpr_is_pos_inf(x) && fmpr_is_neg_inf(y))
        || (fmpr_is_neg_inf(x) && fmpr_is_pos_inf(y)))
    {
        fmpr_nan(z);
        return FMPR_RESULT_EXACT;
    }
    else if (fmpr_is_special(x))
    {
        fmpr_set(z, x);
        return FMPR_RESULT_EXACT;
    }
    else
    {
        fmpr_set(z, y);
        return FMPR_RESULT_EXACT;
    }
}

long
fmpr_add(fmpr_t z, const fmpr_t x, const fmpr_t y, long prec, fmpr_rnd_t rnd)
{
    long shift, xn, yn;
    mp_limb_t xtmp, ytmp;
    mp_ptr xptr, yptr;
    fmpz xv, yv;
    const fmpz * xexp;
    const fmpz * yexp;
    int xsign, ysign;

    if (fmpr_is_special(x) || fmpr_is_special(y))
    {
        return _fmpr_add_special(z, x, y, prec, rnd);
    }

    shift = _fmpz_sub_small(fmpr_expref(y), fmpr_expref(x));

    if (shift >= 0)
    {
        xexp = fmpr_expref(x);
        yexp = fmpr_expref(y);
        xv = *fmpr_manref(x);
        yv = *fmpr_manref(y);
    }
    else
    {
        xexp = fmpr_expref(y);
        yexp = fmpr_expref(x);
        xv = *fmpr_manref(y);
        yv = *fmpr_manref(x);
        shift = -shift;
    }

    READ_SIGN_PTR(xsign, xtmp, xptr, xn, xv)
    READ_SIGN_PTR(ysign, ytmp, yptr, yn, yv)

    if ((xn == 1) && (yn == 1) && (shift < FLINT_BITS))
        return _fmpr_add_small(z, xptr[0], xsign, xexp, yptr[0], ysign, yexp, shift, prec, rnd);
    else
        return _fmpr_add_large(z, xptr, xn, xsign, xexp, yptr, yn, ysign, yexp, shift, prec, rnd);
}


