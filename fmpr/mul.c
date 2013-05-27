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

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "fmpr.h"

#define MUL_STACK_ALLOC 40

/* requires xn >= yn, xman and yman both normalised and odd */
static long
_fmpr_mul_large(fmpz_t zman, fmpz_t zexp,
    mp_srcptr xman, mp_size_t xn, const fmpz_t xexp,
    mp_srcptr yman, mp_size_t yn, const fmpz_t yexp,
    int negative, long prec, fmpr_rnd_t rnd)
{
    long zn, alloc, ret, shift;
    mp_limb_t tmp_stack[MUL_STACK_ALLOC];
    mp_ptr tmp;

    zn = xn + yn;
    alloc = zn;

    if (alloc > MUL_STACK_ALLOC)
        tmp = flint_malloc(alloc * sizeof(mp_limb_t));
    else
        tmp = tmp_stack;

    mpn_mul(tmp, xman, xn, yman, yn);
    zn = zn - (tmp[zn-1] == 0);

    ret = _fmpr_set_round_mpn(&shift, zman, tmp, zn, negative, prec, rnd);
    fmpz_add2_fmpz_si_inline(zexp, xexp, yexp, shift);

    if (alloc > MUL_STACK_ALLOC)
        flint_free(tmp);

    return ret;
}

static void
_fmpr_mul_special(fmpr_t z, const fmpr_t x, const fmpr_t y)
{
    if (fmpr_is_zero(x))
    {
        if (!fmpr_is_special(y) || fmpr_is_zero(y))
            fmpr_zero(z);
        else
            fmpr_nan(z);
        return;
    }

    if (fmpr_is_zero(y))
    {
        if (!fmpr_is_special(x))
            fmpr_zero(z);
        else
            fmpr_nan(z);
        return;
    }

    if ((fmpr_is_inf(x) && (fmpr_is_inf(y) || !fmpr_is_special(y))) ||
        (fmpr_is_inf(y) && !fmpr_is_special(x)))
    {
        if (fmpr_sgn(x) == fmpr_sgn(y))
            fmpr_pos_inf(z);
        else
            fmpr_neg_inf(z);
        return;
    }

    fmpr_nan(z);
}

long
fmpr_mul(fmpr_t z, const fmpr_t x, const fmpr_t y, long prec, fmpr_rnd_t rnd)
{
    __mpz_struct *xmpz, *ympz;
    long xn, yn;
    int negative;
    fmpz u, v;

    if (fmpr_is_special(x) || fmpr_is_special(y))
    {
        _fmpr_mul_special(z, x, y);
        return FMPR_RESULT_EXACT;
    }

    u = *fmpr_manref(x);
    v = *fmpr_manref(y);

    if (!COEFF_IS_MPZ(u))
    {
        if (!COEFF_IS_MPZ(v))
        {
            mp_limb_t hi, lo;
            negative = (u ^ v) < 0;
            u = FLINT_ABS(u);
            v = FLINT_ABS(v);
            umul_ppmm(hi, lo, u, v);

            if (hi == 0)
            {
                /* 1 limb */
                long lead, trail, bc, shift, ret;

                count_leading_zeros(lead, lo);
                bc = FLINT_BITS - lead;

                shift = 0;
                ret = FMPR_RESULT_EXACT;
                if (bc > prec)
                {
                    shift = bc - prec;
                    lo = (lo >> shift) + rounds_up(rnd, negative);
                    count_trailing_zeros(trail, lo);
                    lo >>= trail;
                    shift += trail;
                    ret = trail;
                }

                if (!negative)
                    fmpz_set_ui(fmpr_manref(z), lo);
                else
                    fmpz_neg_ui(fmpr_manref(z), lo);

                fmpz_add2_fmpz_si_inline(fmpr_expref(z), fmpr_expref(x), fmpr_expref(y), shift);

                return ret;
            }
            else
            {
                /* 2 limbs */
                long lead, trail, bc, shift, ret;

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

                fmpz_add2_fmpz_si_inline(fmpr_expref(z), fmpr_expref(x), fmpr_expref(y), shift);
                return ret;
            }
        }
        else
        {
            mp_limb_t t;
            ympz = COEFF_TO_PTR(v);
            yn = ympz->_mp_size;
            negative = (yn < 0) ^ (u < 0);
            t = FLINT_ABS(u);
            yn = FLINT_ABS(yn);
            return _fmpr_mul_large(fmpr_manref(z), fmpr_expref(z),
                ympz->_mp_d, yn, fmpr_expref(y),
                &t, 1, fmpr_expref(x), negative, prec, rnd);
        }
    }
    else
    {
        xmpz = COEFF_TO_PTR(u);
        xn = xmpz->_mp_size;

        if (!COEFF_IS_MPZ(v))
        {
            mp_limb_t t;
            negative = (xn < 0) ^ (v < 0);
            t = FLINT_ABS(v);
            xn = FLINT_ABS(xn);
            return _fmpr_mul_large(fmpr_manref(z), fmpr_expref(z),
                xmpz->_mp_d, xn, fmpr_expref(x),
                &t, 1, fmpr_expref(y), negative, prec, rnd);
        }
        else
        {

            ympz = COEFF_TO_PTR(v);
            yn = ympz->_mp_size;
            negative = (xn ^ yn) < 0;

            xn = FLINT_ABS(xn);
            yn = FLINT_ABS(yn);

            if (xn >= yn)
                return _fmpr_mul_large(fmpr_manref(z), fmpr_expref(z),
                    xmpz->_mp_d, xn, fmpr_expref(x),
                    ympz->_mp_d, yn, fmpr_expref(y), negative, prec, rnd);
            else
                return _fmpr_mul_large(fmpr_manref(z), fmpr_expref(z),
                    ympz->_mp_d, yn, fmpr_expref(y),
                    xmpz->_mp_d, xn, fmpr_expref(x), negative, prec, rnd);
        }
    }
}
