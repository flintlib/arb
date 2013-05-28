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

static long _fmpr_sub_special(fmpr_t z, const fmpr_t x, const fmpr_t y, long prec, fmpr_rnd_t rnd)
{
    if (fmpr_is_zero(x))
    {
        return fmpr_neg_round(z, y, prec, rnd);
    }
    else if (fmpr_is_zero(y))
    {
        return fmpr_set_round(z, x, prec, rnd);
    }
    else if (fmpr_is_nan(x) || fmpr_is_nan(y)
        || (fmpr_is_pos_inf(x) && fmpr_is_pos_inf(y))
        || (fmpr_is_neg_inf(x) && fmpr_is_neg_inf(y)))
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
        fmpr_neg(z, y);
        return FMPR_RESULT_EXACT;
    }
}

long
_fmpr_add_large(fmpr_t z,
        mp_srcptr xman, mp_size_t xn, int xsign, const fmpz_t xexp,
        mp_srcptr yman, mp_size_t yn, int ysign, const fmpz_t yexp,
        long shift, long prec, fmpr_rnd_t rnd);

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


long
fmpr_sub(fmpr_t z, const fmpr_t x, const fmpr_t y, long prec, fmpr_rnd_t rnd)
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
        return _fmpr_sub_special(z, x, y, prec, rnd);
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
    }

    READ_SIGN_PTR(xsign, xtmp, xptr, xn, xv)
    READ_SIGN_PTR(ysign, ytmp, yptr, yn, yv)

    if (shift >= 0)
    {
        ysign = !ysign;
    }
    else
    {
        shift = -shift;
        xsign = !xsign;
    }

    if ((xn == 1) && (yn == 1) && (shift < FLINT_BITS))
        return _fmpr_add_small(z, xptr[0], xsign, xexp, yptr[0], ysign, yexp, shift, prec, rnd);
    else
        return _fmpr_add_large(z, xptr, xn, xsign, xexp, yptr, yn, ysign, yexp, shift, prec, rnd);
}


long
fmpr_sub_ui(fmpr_t z, const fmpr_t x, ulong y, long prec, fmpr_rnd_t rnd)
{
    fmpr_t t; long r;
    fmpr_init(t);
    fmpr_set_ui(t, y);
    r = fmpr_sub(z, x, t, prec, rnd);
    fmpr_clear(t);
    return r;
}

long
fmpr_sub_si(fmpr_t z, const fmpr_t x, long y, long prec, fmpr_rnd_t rnd)
{
    fmpr_t t; long r;
    fmpr_init(t);
    fmpr_set_si(t, y);
    r = fmpr_sub(z, x, t, prec, rnd);
    fmpr_clear(t);
    return r;
}

long
fmpr_sub_fmpz(fmpr_t z, const fmpr_t x, const fmpz_t y, long prec, fmpr_rnd_t rnd)
{
    fmpr_t t; long r;
    fmpr_init(t);
    fmpr_set_fmpz(t, y);
    r = fmpr_sub(z, x, t, prec, rnd);
    fmpr_clear(t);
    return r;
}
