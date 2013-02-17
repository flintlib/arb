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

#include "fmprb.h"

/*
We speed up the radius operations by working with mantissas aligned to
FMPRB_RAD_PREC bits (possibly one less bit after multiplying).

The error is computed as x*b + y*a + a*b + r where x, b, y, a, a, b, r are
floating-point numbers of the form m * 2^e and m has
exactly FMPRB_RAD_PREC bits, i.e. m <= 2^FMPRB_RAD_PREC - 1.

The mantissas of the products x*b, y*a, a*b are computed as
m3 = floor((m1 * m2) / 2^FMPRB_RAD_PREC) + 1.
One can then verify that m3 <= 2^FMPRB_RAD_PREC - 1.

In the additions, we do not normalise the output mantissa. The result is
largest when exponents are the same: then we have precisely m3 = m1 + m2.

Thus the final mantissa is <= 4 * (2^FMPRB_RAD_PREC - 1).
If FMPRB_RAD_PREC <= FLINT_BITS - 2, this precisely fits in a limb
without overflow.
*/

#define _RAD_MUL(a, ae, b, be) \
    do { \
        mp_limb_t hi, lo; \
        umul_ppmm(hi, lo, a, b); \
        a = ((hi << (FLINT_BITS - FMPRB_RAD_PREC)) | (lo >> FMPRB_RAD_PREC)) + 1; \
        fmpz_add_inline(ae, ae, be); \
    } while (0); \

static __inline__ void
_rad_bound(mp_limb_t * m, fmpz_t exp, const fmpr_t t)
{
    long e;
    *m = fmpz_abs_ubound_ui_2exp(&e, fmpr_manref(t), FMPRB_RAD_PREC);
    fmpz_add_si_inline(exp, fmpr_expref(t), e);
}

static __inline__ mp_limb_t
_rad_add(mp_limb_t a, fmpz_t aexp, mp_limb_t b, fmpz_t bexp)
{
    long shift;

    shift = _fmpz_sub_small(aexp, bexp);

    if (shift == 0)
    {
        return a + b;
    }
    else if (shift > 0)
    {
        if (shift <= FMPRB_RAD_PREC)
            return a + (b >> shift) + 1;
        else
            return a + 1;
    }
    else
    {
        fmpz_swap(aexp, bexp);

        if ((-shift) <= FMPRB_RAD_PREC)
            return b + (a >> (-shift)) + 1;
        else
            return b + 1;
    }
}

void _fmprb_mul_main(fmpr_t z, fmpr_t c,
    const fmpr_t x, const fmpr_t a,
    const fmpr_t y, const fmpr_t b, long prec)
{
    mp_limb_t xm, am, ym, bm;
    fmpz_t xe, ae, ye, be;
    long r, shift;

    fmpz_init(xe);
    fmpz_init(ae);
    fmpz_init(ye);
    fmpz_init(be);

    _rad_bound(&xm, xe, x);
    _rad_bound(&ym, ye, y);
    _rad_bound(&am, ae, a);
    _rad_bound(&bm, be, b);

    _RAD_MUL(xm, xe, bm, be);
    _RAD_MUL(ym, ye, am, ae);
    xm = _rad_add(xm, xe, ym, ye);

    _RAD_MUL(am, ae, bm, be);
    xm = _rad_add(xm, xe, am, ae);

    r = fmpr_mul(z, x, y, prec, FMPR_RND_DOWN);

    if (r != FMPR_RESULT_EXACT)
    {
        am = 1UL << FMPRB_RAD_PREC;
        fmpz_add_si_inline(ae, fmpr_expref(z), -r - 2*FMPRB_RAD_PREC);
        xm = _rad_add(xm, xe, am, ae);
    }

    shift = FMPRB_RAD_PREC;

    /* make the radius mantissa odd and small */
    xm += !(xm & 1);
    while (xm >= (1UL << FMPRB_RAD_PREC))
    {
        xm = (xm >> 1) + 1;
        xm += !(xm & 1);
        shift++;
    }

    fmpz_set_ui(fmpr_manref(c), xm);
    fmpz_add_si_inline(fmpr_expref(c), xe, shift);

    fmpz_clear(xe);
    fmpz_clear(ae);
    fmpz_clear(ye);
    fmpz_clear(be);
}

void _fmprb_mul_fmpr_main(fmpr_t z, fmpr_t c,
    const fmpr_t x, const fmpr_t a,
    const fmpr_t y, long prec)
{
    mp_limb_t ym, am;
    fmpz_t ye, ae;
    long r, shift;

    fmpz_init(ye);
    fmpz_init(ae);

    _rad_bound(&ym, ye, y);
    _rad_bound(&am, ae, a);
    _RAD_MUL(ym, ye, am, ae);

    r = fmpr_mul(z, x, y, prec, FMPR_RND_DOWN);

    if (r != FMPR_RESULT_EXACT)
    {
        am = 1UL << FMPRB_RAD_PREC;
        fmpz_add_si_inline(ae, fmpr_expref(z), -r - 2*FMPRB_RAD_PREC);
        ym = _rad_add(ym, ye, am, ae);
    }

    shift = FMPRB_RAD_PREC;

    /* make the radius mantissa odd and small */
    ym += !(ym & 1);
    while (ym >= (1UL << FMPRB_RAD_PREC))
    {
        ym = (ym >> 1) + 1;
        ym += !(ym & 1);
        shift++;
    }

    fmpz_set_ui(fmpr_manref(c), ym);
    fmpz_add_si_inline(fmpr_expref(c), ye, shift);

    fmpz_clear(ye);
    fmpz_clear(ae);
}

void
fmprb_mul_fmpr(fmprb_t z, const fmprb_t x, const fmpr_t y, long prec)
{
    if (fmprb_is_exact(x))
    {
        long r;
        r = fmpr_mul(fmprb_midref(z), fmprb_midref(x), y, prec, FMPR_RND_DOWN);
        fmpr_set_error_result(fmprb_radref(z), fmprb_midref(z), r);
    }
    else
    {
        if (fmpr_is_special(fmprb_midref(x)) ||
            fmpr_is_special(fmprb_radref(x)) || fmpr_is_special(y))
        {
            fmprb_mul_fmpr_naive(z, x, y, prec);
        }
        else
        {
            _fmprb_mul_fmpr_main(fmprb_midref(z), fmprb_radref(z),
                fmprb_midref(x), fmprb_radref(x), y, prec);
        }
    }
}

void
fmprb_mul(fmprb_t z, const fmprb_t x, const fmprb_t y, long prec)
{
    if (fmprb_is_exact(x))
    {
        fmprb_mul_fmpr(z, y, fmprb_midref(x), prec);
    }
    else if (fmprb_is_exact(y))
    {
        fmprb_mul_fmpr(z, x, fmprb_midref(y), prec);
    }
    else
    {
        if (fmpr_is_special(fmprb_midref(x)) || 
                fmpr_is_special(fmprb_radref(x)) ||
                fmpr_is_special(fmprb_midref(y)) ||
                fmpr_is_special(fmprb_radref(y)))
        {
            fmprb_mul_main_naive(z, x, y, prec);
        }
        else
        {
            _fmprb_mul_main(fmprb_midref(z), fmprb_radref(z),
                fmprb_midref(x), fmprb_radref(x),
                fmprb_midref(y), fmprb_radref(y), prec);
        }
    }
}

void
fmprb_mul_ui(fmprb_t z, const fmprb_t x, ulong y, long prec)
{
    fmpr_t t;
    fmpr_init(t);
    fmpr_set_ui(t, y);
    fmprb_mul_fmpr(z, x, t, prec);
    fmpr_clear(t);
}

void
fmprb_mul_si(fmprb_t z, const fmprb_t x, long y, long prec)
{
    fmpr_t t;
    fmpr_init(t);
    fmpr_set_si(t, y);
    fmprb_mul_fmpr(z, x, t, prec);
    fmpr_clear(t);
}

void
fmprb_mul_fmpz(fmprb_t z, const fmprb_t x, const fmpz_t y, long prec)
{
    fmpr_t t;
    fmpr_init(t);
    fmpr_set_fmpz(t, y);
    fmprb_mul_fmpr(z, x, t, prec);
    fmpr_clear(t);
}

