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

    Copyright (C) 2013 Fredrik Johansson

******************************************************************************/

#ifndef FMPZ_EXTRAS_H
#define FMPZ_EXTRAS_H

#include <limits.h>
#include "flint.h"
#include "fmpz.h"

#ifdef __cplusplus
extern "C" {
#endif

static __inline__ void
fmpz_add_inline(fmpz_t z, const fmpz_t x, const fmpz_t y)
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
fmpz_add_si(fmpz_t z, const fmpz_t x, long y)
{
    if (y >= 0)
        fmpz_add_ui(z, x, y);
    else
        fmpz_sub_ui(z, x, -y);
}

static __inline__ void
fmpz_sub_si(fmpz_t z, const fmpz_t x, long y)
{
    if (y >= 0)
        fmpz_sub_ui(z, x, y);
    else
        fmpz_add_ui(z, x, -y);
}

static __inline__ void
fmpz_add_si_inline(fmpz_t z, const fmpz_t x, long y)
{
    fmpz f;

    f = *x;

    if (!COEFF_IS_MPZ(f) && (COEFF_MIN <= y && y <= COEFF_MAX))
        fmpz_set_si(z, f + y);
    else
        fmpz_add_si(z, x, y);
}

static __inline__ void
fmpz_sub_si_inline(fmpz_t z, const fmpz_t x, long y)
{
    fmpz f;

    f = *x;

    if (!COEFF_IS_MPZ(f) && (COEFF_MIN <= y && y <= COEFF_MAX))
        fmpz_set_si(z, f - y);
    else
        fmpz_sub_si(z, x, y);
}

static __inline__ void
fmpz_add_ui_inline(fmpz_t z, const fmpz_t x, ulong y)
{
    fmpz f = *x;

    if (!COEFF_IS_MPZ(f) && y <= COEFF_MAX)
        fmpz_set_si(z, f + y);
    else
        fmpz_add_ui(z, x, y);
}

static __inline__ void
fmpz_add2_fmpz_si_inline(fmpz_t z, const fmpz_t x, const fmpz_t y, long c)
{
    fmpz f, g, h;

    f = *x;
    g = *y;

    if (!COEFF_IS_MPZ(f) && !COEFF_IS_MPZ(g))
    {
        h = f + g;

        if ((COEFF_MIN <= h && h <= COEFF_MAX) &&
            (COEFF_MIN <= c && c <= COEFF_MAX))
        {
            fmpz_set_si(z, h + c);
            return;
        }
    }

    fmpz_add(z, x, y);
    fmpz_add_si(z, z, c);
}

/* sets z = +/- ({src, n} >> shift) where 0 <= shift < FLINT_BITS
   and the top limb of src is nonzero */
/* TODO: optimize for result = 1 limb */

static __inline__ void
fmpz_set_mpn_rshift(fmpz_t z, mp_srcptr src, mp_size_t n, unsigned int shift, int negative)
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
fmpz_set_mpn_large(fmpz_t z, mp_srcptr src, mp_size_t n, int negative)
{
    __mpz_struct * zz;
    zz = _fmpz_promote(z);

    if (zz->_mp_alloc < n)
        mpz_realloc2(zz, n * FLINT_BITS);

    flint_mpn_copyi(zz->_mp_d, src, n);
    zz->_mp_size = negative ? -n : n;
}

/* round away from zero */
static __inline__ void
fmpz_adiv_q_2exp(fmpz_t z, const fmpz_t x, mp_bitcnt_t exp)
{
    int sign = fmpz_sgn(x);

    if (sign > 0)
        fmpz_cdiv_q_2exp(z, x, exp);
    else
        fmpz_fdiv_q_2exp(z, x, exp);
}

/* sets z = x + y*2^shift */
static __inline__ void fmpz_add_mul2exp(fmpz_t z, const fmpz_t x, const fmpz_t y, ulong shift)
{
    fmpz_t t;
    fmpz_init(t);
    fmpz_mul_2exp(t, y, shift);
    fmpz_add(z, x, t);
    fmpz_clear(t);
}

static __inline__ void fmpz_sub_mul2exp(fmpz_t z, const fmpz_t x, const fmpz_t y, ulong shift)
{
    fmpz_t t;
    fmpz_init(t);
    fmpz_mul_2exp(t, y, shift);
    fmpz_sub(z, x, t);
    fmpz_clear(t);
}

long _fmpz_sub_small_large(const fmpz_t x, const fmpz_t y);

static __inline__ long
_fmpz_sub_small(const fmpz_t x, const fmpz_t y)
{
    if (!COEFF_IS_MPZ(*x) && !COEFF_IS_MPZ(*y))
    {
        return (*x) - (*y);
    }
    else
    {
        return _fmpz_sub_small_large(x, y);
    }
}

static __inline__ mp_size_t
_fmpz_size(const fmpz_t f)
{
    fmpz d = *f;

    if (!COEFF_IS_MPZ(d))
        return 1;
    else
        return mpz_size(COEFF_TO_PTR(d));
}

static __inline__ void
fmpz_ui_pow_ui(fmpz_t x, ulong b, ulong e)
{
    if (e <= 1)
    {
        fmpz_set_ui(x, e == 0 ? 1UL : b);
    }
    else
    {
        fmpz_set_ui(x, b);
        fmpz_pow_ui(x, x, e);
    }
}

static __inline__ void
fmpz_max(fmpz_t z, const fmpz_t x, const fmpz_t y)
{
    if (fmpz_cmp(x, y) >= 0)
        fmpz_set(z, x);
    else
        fmpz_set(z, y);
}

static __inline__ void
fmpz_min(fmpz_t z, const fmpz_t x, const fmpz_t y)
{
    if (fmpz_cmp(x, y) < 0)
        fmpz_set(z, x);
    else
        fmpz_set(z, y);
}

#define FMPZ_GET_MPN_READONLY(zsign, zn, zptr, ztmp, zv) \
    if (!COEFF_IS_MPZ(zv)) \
    { \
        (zsign) = (zv) < 0; \
        (ztmp) = FLINT_ABS(zv); \
        (zptr) = &(ztmp); \
        (zn) = 1; \
    } \
    else \
    { \
        __mpz_struct * ___zz = COEFF_TO_PTR(zv); \
        (zptr) = ___zz->_mp_d; \
        (zn) = ___zz->_mp_size; \
        (zsign) = (zn) < 0; \
        (zn) = FLINT_ABS(zn); \
    }

void fmpz_lshift_mpn(fmpz_t z, mp_srcptr d, mp_size_t dn, int sgnbit, mp_bitcnt_t shift);

#ifdef __cplusplus
}
#endif

#endif

