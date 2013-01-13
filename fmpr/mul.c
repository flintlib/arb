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

static __inline__ int
rounds_up(fmpr_rnd_t rnd, int negative)
{
    if (rnd == FMPR_RND_DOWN) return 0;
    if (rnd == FMPR_RND_UP) return 1;
    if (rnd == FMPR_RND_FLOOR) return negative;
    return !negative;
}

static __inline__ void
_fmpz_add(fmpz_t z, const fmpz_t x, const fmpz_t y)
{
    fmpz f, g;

    f = *x;
    g = *y;

    if (!COEFF_IS_MPZ(f) && !COEFF_IS_MPZ(g))
        fmpz_set_si(z, f + g);
    else
        fmpz_add(z, x, y);
}

static __inline__ void
_fmpz_add2_fmpz_si(fmpz_t z, const fmpz_t x, const fmpz_t y, long c)
{
    fmpz f, g, h;

    f = *x;
    g = *y;

    if (!COEFF_IS_MPZ(f) && !COEFF_IS_MPZ(g))
    {
        h = f + g;

        if (!COEFF_IS_MPZ(h) && !COEFF_IS_MPZ(c))
        {
            fmpz_set_si(z, h + c);
            return;
        }
    }

    fmpz_add(z, x, y);
    if (c >= 0)
        fmpz_add_ui(z, z, c);
    else
        fmpz_sub_ui(z, z, -c);
}

/* sets z = +/- ({src, n} >> shift) where 0 <= shift < FLINT_BITS
   and the top limb of src is nonzero */
/* TODO: optimize for result = 1 limb */

static __inline__ void
_fmpz_set_mpn_rshift(fmpz_t z, mp_ptr src, mp_size_t n, unsigned int shift, int negative)
{
    __mpz_struct * zptr;
    zptr = _fmpz_promote(z);

    if (zptr->_mp_alloc < n)
        mpz_realloc2(zptr, n * FLINT_BITS);

    if (shift == 0)
    {
        flint_mpn_copyi(zptr->_mp_d, src, n);
    }
    else
    {
        mpn_rshift(zptr->_mp_d, src, n, shift);
        while (zptr->_mp_d[n - 1] == 0) /* todo: can only do one iter? */
            n--;
    }

    zptr->_mp_size = negative ? -n : n;
    _fmpz_demote_val(z);
}

/* sets z = +/- {src, n} where n >= 2 and the top limb of src is nonzero */
static __inline__ void
_fmpz_set_mpn_large(fmpz_t z, mp_ptr src, mp_size_t n, int negative)
{
    __mpz_struct * zz;
    zz = _fmpz_promote(z);

    if (zz->_mp_alloc < n)
        mpz_realloc2(zz, n * FLINT_BITS);

    flint_mpn_copyi(zz->_mp_d, src, n);
    zz->_mp_size = negative ? -n : n;
}

#define MUL_STACK_ALLOC 40

/* requires xn >= yn, xman and yman both normalised and odd */
static long
_fmpr_mul_large(fmpz_t zman, fmpz_t zexp,
    mp_srcptr xman, mp_size_t xn, const fmpz_t xexp,
    mp_srcptr yman, mp_size_t yn, const fmpz_t yexp,
    int negative, long prec, fmpr_rnd_t rnd)
{
    long zn, alloc, lead, bc, ret;
    mp_limb_t tmp_stack[MUL_STACK_ALLOC];
    mp_ptr tmp;

    zn = xn + yn;
    alloc = zn + 1;  /* may need one extra temp limb for rounding carry */

    if (alloc > MUL_STACK_ALLOC)
        tmp = flint_malloc(alloc * sizeof(mp_limb_t));
    else
        tmp = tmp_stack;

    mpn_mul(tmp, xman, xn, yman, yn);

    zn = zn - (tmp[zn-1] == 0);
    count_leading_zeros(lead, tmp[zn-1]);
    bc = zn * FLINT_BITS - lead;

    if (bc <= prec)
    {
        if (zn == 1)
        {
            if (!negative)
                fmpz_set_ui(zman, tmp[0]);
            else
                fmpz_neg_ui(zman, tmp[0]);
        }
        else
        {
            _fmpz_set_mpn_large(zman, tmp, zn, negative);
        }

        _fmpz_add(zexp, xexp, yexp);
        ret = FMPR_RESULT_EXACT;
    }
    else
    {
        long cut_limb, cut_bit, shift, trail_limbs, trail_bits;
        mp_limb_t cy;
        mp_ptr res;

        shift = bc - prec;
        cut_limb = shift / FLINT_BITS;
        cut_bit = shift % FLINT_BITS;

        /* cut off at the rounding bit */
        tmp[cut_limb] &= ~(((mp_limb_t) 1UL << cut_bit) - 1UL);
        if (rounds_up(rnd, negative))
        {
            cy = mpn_add_1(tmp + cut_limb,
                tmp + cut_limb, zn - cut_limb, (1UL << cut_bit));
            tmp[zn] = cy;
            zn += cy;
        }

        /* remove trailing zero limbs */
        trail_limbs = 0;
        while (tmp[cut_limb + trail_limbs] == 0)
            trail_limbs++;

        res = tmp + (cut_limb + trail_limbs);
        zn -= (cut_limb + trail_limbs);
        count_trailing_zeros(trail_bits, res[0]);
        ret = (trail_limbs * FLINT_BITS) + trail_bits - cut_bit;
        shift = (cut_limb + trail_limbs) * FLINT_BITS + trail_bits;

        _fmpz_set_mpn_rshift(zman, res, zn, trail_bits, negative);
        _fmpz_add2_fmpz_si(zexp, xexp, yexp, shift);
    }

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

                _fmpz_add2_fmpz_si(fmpr_expref(z), fmpr_expref(x), fmpr_expref(y), shift);

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

                _fmpz_add2_fmpz_si(fmpr_expref(z), fmpr_expref(x), fmpr_expref(y), shift);
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
