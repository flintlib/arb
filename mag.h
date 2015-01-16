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

#ifndef MAG_H
#define MAG_H

#ifdef MAG_INLINES_C
#define MAG_INLINE
#else
#define MAG_INLINE static __inline__
#endif

#include <math.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_extras.h"

#ifdef __cplusplus
extern "C" {
#endif

#define LIMB_ONE ((mp_limb_t) 1)
#define LIMB_ONES (-(mp_limb_t) 1)
#define LIMB_TOP (((mp_limb_t) 1) << (FLINT_BITS - 1))
#define MASK_LIMB(n, c) ((n) & (LIMB_ONES << (c)))

#define MAG_MAX_LAGOM_EXP (COEFF_MAX / 4)
#define MAG_MIN_LAGOM_EXP (-MAG_MAX_LAGOM_EXP)

#define ADD2_FAST_MAX (COEFF_MAX / 4)
#define ADD2_FAST_MIN (-ADD2_FAST_MAX)

/* TODO: rename these and move them to fmpz_extras */

static __inline__ void
_fmpz_set_fast(fmpz_t f, const fmpz_t g)
{
    if (!COEFF_IS_MPZ(*f) && !COEFF_IS_MPZ(*g))
        *f = *g;
    else
        fmpz_set(f, g);
}

static __inline__ void
_fmpz_add_fast(fmpz_t z, const fmpz_t x, long c)
{
    fmpz ze, xe;

    ze = *z;
    xe = *x;

    if (!COEFF_IS_MPZ(ze) && (xe > ADD2_FAST_MIN && xe < ADD2_FAST_MAX))
        *z = xe + c;
    else
        fmpz_add_si(z, x, c);
}

static __inline__ void
_fmpz_add2_fast(fmpz_t z, const fmpz_t x, const fmpz_t y, long c)
{
    fmpz ze, xe, ye;

    ze = *z;
    xe = *x;
    ye = *y;

    if (!COEFF_IS_MPZ(ze) && (xe > ADD2_FAST_MIN && xe < ADD2_FAST_MAX) &&
                             (ye > ADD2_FAST_MIN && ye < ADD2_FAST_MAX))
    {
        *z = xe + ye + c;
    }
    else
    {
        fmpz_add(z, x, y);
        fmpz_add_si(z, z, c);
    }
}

static __inline__ void
_fmpz_sub2_fast(fmpz_t z, const fmpz_t x, const fmpz_t y, long c)
{
    fmpz ze, xe, ye;

    ze = *z;
    xe = *x;
    ye = *y;

    if (!COEFF_IS_MPZ(ze) && (xe > ADD2_FAST_MIN && xe < ADD2_FAST_MAX) &&
                             (ye > ADD2_FAST_MIN && ye < ADD2_FAST_MAX))
    {
        *z = xe - ye + c;
    }
    else
    {
        fmpz_sub(z, x, y);
        fmpz_add_si(z, z, c);
    }
}


#define MAG_EXP_POS_INF (COEFF_MIN+1)

/* Finite and with lagom big exponents. */
#define MAG_IS_LAGOM(x) (MAG_EXP(x) >= MAG_MIN_LAGOM_EXP && \
                         MAG_EXP(x) <= MAG_MAX_LAGOM_EXP)

#define MAG_EXPREF(x) (&(x)->exp)
#define MAG_EXP(x) ((x)->exp)
#define MAG_MAN(x) ((x)->man)

#define MAG_BITS 30

#define MAG_ONE_HALF (1UL << (MAG_BITS - 1))

static __inline__ mp_limb_t
__mag_fixmul32(mp_limb_t x, mp_limb_t y)
{
    mp_limb_t u, v;
    umul_ppmm(u, v, x, y);
    return (u << (32 - MAG_BITS)) | (v >> MAG_BITS);
}

#if FLINT_BITS == 64
#define MAG_FIXMUL(x, y) (((x) * (y)) >> MAG_BITS)
#else
#define MAG_FIXMUL(x, y) __mag_fixmul32((x), (y))
#endif

#define MAG_CHECK_BITS(rr) \
    if (MAG_MAN(rr) != 0 && FLINT_BIT_COUNT(MAG_MAN(rr)) != MAG_BITS) \
    { \
        printf("FAIL: wrong number of bits in mantissa!\n"); \
        abort(); \
    }

/* Note: assumes mantissa either has the right number of bits, or
   one more bit (but in that case must not be all ones, as that
   would round up to one extra bit again, requiring a second
   correction). */

#define MAG_ADJUST_ONE_TOO_LARGE(x) \
    do { \
        mp_limb_t __t = MAG_MAN(x) >> MAG_BITS; \
        MAG_MAN(x) = (MAG_MAN(x) >> __t) + __t; \
        if (__t) \
            fmpz_add_ui(MAG_EXPREF(x), MAG_EXPREF(x), __t); \
    } while (0)

#define MAG_FAST_ADJUST_ONE_TOO_LARGE(x) \
    do { \
        mp_limb_t __t = MAG_MAN(x) >> MAG_BITS; \
        MAG_MAN(x) = (MAG_MAN(x) >> __t) + __t; \
        MAG_EXP(x) += __t; \
    } while (0)

#define MAG_ADJUST_ONE_TOO_SMALL(x) \
    do { \
        mp_limb_t __t = !(MAG_MAN(x) >> (MAG_BITS - 1)); \
        MAG_MAN(x) = (MAG_MAN(x) << __t); \
        if (__t) \
            fmpz_sub_ui(MAG_EXPREF(x), MAG_EXPREF(x), __t); \
    } while (0)

#define MAG_FAST_ADJUST_ONE_TOO_SMALL(x) \
    do { \
        mp_limb_t __t = !(MAG_MAN(x) >> (MAG_BITS - 1)); \
        MAG_MAN(x) = (MAG_MAN(x) << __t); \
        MAG_EXP(x) -= __t; \
    } while (0)


typedef struct
{
    fmpz exp;
    mp_limb_t man;
}
mag_struct;

typedef mag_struct mag_t[1];
typedef mag_struct * mag_ptr;
typedef const mag_struct * mag_srcptr;

MAG_INLINE void
mag_init(mag_t x)
{
    fmpz_init(MAG_EXPREF(x));
    MAG_MAN(x) = 0;
}

MAG_INLINE void
mag_init_set(mag_t x, const mag_t y)
{
    fmpz_init_set(MAG_EXPREF(x), MAG_EXPREF(y));
    MAG_MAN(x) = MAG_MAN(y);
}

MAG_INLINE void
mag_clear(mag_t x)
{
    fmpz_clear(MAG_EXPREF(x));
}

MAG_INLINE void
mag_swap(mag_t x, mag_t y)
{
    mag_struct t = *x;
    *x = *y;
    *y = t;
}

MAG_INLINE void
mag_set(mag_t x, const mag_t y)
{
    _fmpz_set_fast(MAG_EXPREF(x), MAG_EXPREF(y));
    x->man = y->man;
}

MAG_INLINE void
mag_zero(mag_t x)
{
    fmpz_zero(MAG_EXPREF(x));
    MAG_MAN(x) = 0;
}

MAG_INLINE void
mag_one(mag_t x)
{
    fmpz_one(MAG_EXPREF(x));
    MAG_MAN(x) = MAG_ONE_HALF;
}

MAG_INLINE int
mag_is_special(const mag_t x)
{
    return MAG_MAN(x) == 0;
}

MAG_INLINE int
mag_is_zero(const mag_t x)
{
    return (MAG_MAN(x) == 0) && (MAG_EXP(x) == 0);
}

MAG_INLINE void
mag_inf(mag_t x)
{
    fmpz_clear(MAG_EXPREF(x));
    MAG_EXP(x) = MAG_EXP_POS_INF;
    MAG_MAN(x) = 0;
}

MAG_INLINE int
mag_is_inf(const mag_t x)
{
    return (MAG_MAN(x) == 0) && (MAG_EXP(x) != 0);
}

MAG_INLINE int
mag_is_finite(const mag_t x)
{
    return !mag_is_inf(x);
}

MAG_INLINE int
mag_equal(const mag_t x, const mag_t y)
{
    return (MAG_MAN(x) == MAG_MAN(y))
        && fmpz_equal(MAG_EXPREF(x), MAG_EXPREF(y));
}

/* general versions */

void mag_mul(mag_t z, const mag_t x, const mag_t y);

void mag_mul_lower(mag_t z, const mag_t x, const mag_t y);

void mag_addmul(mag_t z, const mag_t x, const mag_t y);

void mag_add_2exp_fmpz(mag_t z, const mag_t x, const fmpz_t e);

void mag_add(mag_t z, const mag_t x, const mag_t y);

void mag_add_lower(mag_t z, const mag_t x, const mag_t y);

void mag_div(mag_t z, const mag_t x, const mag_t y);

MAG_INLINE void
mag_mul_2exp_si(mag_t z, const mag_t x, long y)
{
    if (mag_is_special(x))
    {
        mag_set(z, x);
    }
    else
    {
        if (y >= ADD2_FAST_MIN && y <= ADD2_FAST_MAX)
            _fmpz_add_fast(MAG_EXPREF(z), MAG_EXPREF(x), y);
        else
            fmpz_add_si(MAG_EXPREF(z), MAG_EXPREF(x), y);
        MAG_MAN(z) = MAG_MAN(x);
    }
}

MAG_INLINE void
mag_mul_2exp_fmpz(mag_t z, const mag_t x, const fmpz_t y)
{
    if (mag_is_special(x))
    {
        mag_set(z, x);
    }
    else
    {
        _fmpz_add2_fast(MAG_EXPREF(z), MAG_EXPREF(x), y, 0);
        MAG_MAN(z) = MAG_MAN(x);
    }
}

void mag_sub_lower(mag_t z, const mag_t x, const mag_t y);

/* Fast versions (no infs/nans, small exponents). Note that this
   applies to outputs too! */

MAG_INLINE void
mag_fast_init_set(mag_t x, const mag_t y)
{
    MAG_EXP(x) = MAG_EXP(y);
    MAG_MAN(x) = MAG_MAN(y);
}

MAG_INLINE void
mag_fast_zero(mag_t x)
{
    MAG_EXP(x) = 0;
    MAG_MAN(x) = 0;
}

MAG_INLINE int
mag_fast_is_zero(const mag_t x)
{
    return MAG_MAN(x) == 0;
}

MAG_INLINE void
mag_fast_mul(mag_t z, const mag_t x, const mag_t y)
{
    if (MAG_MAN(x) == 0 || MAG_MAN(y) == 0)
    {
        mag_fast_zero(z);
    }
    else
    {
        MAG_MAN(z) = MAG_FIXMUL(MAG_MAN(x), MAG_MAN(y)) + LIMB_ONE;
        MAG_EXP(z) = MAG_EXP(x) + MAG_EXP(y);
        MAG_FAST_ADJUST_ONE_TOO_SMALL(z);
    }
}

MAG_INLINE void
mag_fast_mul_2exp_si(mag_t z, const mag_t x, long y)
{
    if (MAG_MAN(x) == 0)
    {
        mag_fast_zero(z);
    }
    else
    {
        MAG_MAN(z) = MAG_MAN(x);
        MAG_EXP(z) = MAG_EXP(x) + y;
    }
}

MAG_INLINE void
mag_fast_addmul(mag_t z, const mag_t x, const mag_t y)
{
    if (MAG_MAN(z) == 0)
    {
        mag_fast_mul(z, x, y);
    }
    else if (MAG_MAN(x) == 0 || MAG_MAN(y) == 0)
    {
        return;
    }
    else
    {
        long shift, e;

        /* x*y < 2^e */
        e = MAG_EXP(x) + MAG_EXP(y);
        shift = MAG_EXP(z) - e;

        if (shift >= 0)
        {
            if (shift >= MAG_BITS)
                MAG_MAN(z)++;
            else
                MAG_MAN(z) = MAG_MAN(z) + (MAG_FIXMUL(MAG_MAN(x), MAG_MAN(y)) >> shift) + 1;
        }
        else
        {
            shift = -shift;
            MAG_EXP(z) = e;

            if (shift >= MAG_BITS)
                MAG_MAN(z) = MAG_FIXMUL(MAG_MAN(x), MAG_MAN(y)) + 2;
            else
                MAG_MAN(z) = MAG_FIXMUL(MAG_MAN(x), MAG_MAN(y)) + (MAG_MAN(z) >> shift) + 2;

            MAG_FAST_ADJUST_ONE_TOO_SMALL(z);
        }

        MAG_FAST_ADJUST_ONE_TOO_LARGE(z);
    }
}

MAG_INLINE void
mag_fast_add_2exp_si(mag_t z, const mag_t x, long e)
{
    /* Must be zero */
    if (mag_is_special(x))
    {
        MAG_MAN(z) = MAG_ONE_HALF;
        MAG_EXP(z) = e + 1;
    }
    else
    {
        long shift;
        shift = MAG_EXP(x) - e;

        if (shift > 0)
        {
            MAG_EXP(z) = MAG_EXP(x);

            if (shift >= MAG_BITS)
                MAG_MAN(z) = MAG_MAN(x) + LIMB_ONE;
            else
                MAG_MAN(z) = MAG_MAN(x) + (LIMB_ONE << (MAG_BITS - shift));
        }
        else
        {
            shift = -shift;

            MAG_EXP(z) = e + 1;

            if (shift >= MAG_BITS)
                MAG_MAN(z) = MAG_ONE_HALF + LIMB_ONE;
            else
                MAG_MAN(z) = MAG_ONE_HALF + (MAG_MAN(x) >> (shift + 1)) + LIMB_ONE;
        }

        MAG_FAST_ADJUST_ONE_TOO_LARGE(z);
    }
}

void mag_set_d_2exp_fmpz(mag_t z, double c, const fmpz_t exp);
void mag_set_fmpz_2exp_fmpz(mag_t z, const fmpz_t man, const fmpz_t exp);

#include "fmpr.h"

void mag_set_fmpr(mag_t x, const fmpr_t y);

void mag_get_fmpr(fmpr_t x, const mag_t r);

void mag_randtest_special(mag_t x, flint_rand_t state, long expbits);

void mag_randtest(mag_t x, flint_rand_t state, long expbits);

void mag_print(const mag_t x);

void mag_printd(const mag_t x, long d);

void mag_get_fmpq(fmpq_t y, const mag_t x);

int mag_cmp(const mag_t x, const mag_t y);

int mag_cmp_2exp_si(const mag_t x, long e);

MAG_INLINE void
mag_min(mag_t z, const mag_t x, const mag_t y)
{
    if (mag_cmp(x, y) <= 0)
        mag_set(z, x);
    else
        mag_set(z, y);
}

MAG_INLINE void
mag_max(mag_t z, const mag_t x, const mag_t y)
{
    if (mag_cmp(x, y) >= 0)
        mag_set(z, x);
    else
        mag_set(z, y);
}

MAG_INLINE mag_ptr
_mag_vec_init(long n)
{
    long i;
    mag_ptr v = (mag_ptr) flint_malloc(sizeof(mag_struct) * n);

    for (i = 0; i < n; i++)
        mag_init(v + i);

    return v;
}

MAG_INLINE void
_mag_vec_clear(mag_ptr v, long n)
{
    long i;
    for (i = 0; i < n; i++)
        mag_clear(v + i);
    flint_free(v);
}

MAG_INLINE void mag_set_d(mag_t z, double x)
{
    fmpz_t e;
    fmpz_init(e);
    mag_set_d_2exp_fmpz(z, x, e);
    fmpz_clear(e);
}

/* TODO: test/document */
double mag_get_d(const mag_t z);

/* TODO: document */
double mag_d_log_upper_bound(double x);
double mag_d_log_lower_bound(double x);

void mag_log1p(mag_t z, const mag_t x);

void mag_log_ui(mag_t t, ulong n);

void mag_exp_maglim(mag_t y, const mag_t x, long maglim);

MAG_INLINE void
mag_exp(mag_t y, const mag_t x)
{
    mag_exp_maglim(y, x, 128);
}

void mag_expm1(mag_t y, const mag_t x);
void mag_exp_tail(mag_t z, const mag_t x, ulong N);

void mag_pow_ui(mag_t z, const mag_t x, ulong e);
void mag_pow_ui_lower(mag_t z, const mag_t x, ulong e);
void mag_pow_fmpz(mag_t z, const mag_t x, const fmpz_t e);

void mag_fac_ui(mag_t z, ulong n);
void mag_rfac_ui(mag_t z, ulong n);

/* TODO: test */
void mag_bernoulli_div_fac_ui(mag_t z, ulong n);

/* TODO: test */
void mag_set_fmpz_2exp_fmpz_lower(mag_t z, const fmpz_t man, const fmpz_t exp);

void mag_sqrt(mag_t y, const mag_t x);
void mag_rsqrt(mag_t y, const mag_t x);

void mag_hypot(mag_t z, const mag_t x, const mag_t y);

void mag_binpow_uiui(mag_t b, ulong m, ulong n);

void mag_polylog_tail(mag_t u, const mag_t z, long sigma, ulong d, ulong N);

void mag_set_ui(mag_t z, ulong x);
void mag_set_ui_lower(mag_t z, ulong x);

/* TODO: test functions below */
MAG_INLINE void
mag_set_ui_2exp_si(mag_t z, ulong v, long e)
{
    mag_set_ui(z, v);
    mag_mul_2exp_si(z, z, e);
}

MAG_INLINE void
mag_set_fmpz(mag_t z, const fmpz_t x)
{
    fmpz_t exp;
    *exp = 0;
    mag_set_fmpz_2exp_fmpz(z, x, exp);
}

MAG_INLINE void
mag_set_fmpz_lower(mag_t z, const fmpz_t x)
{
    fmpz_t exp;
    *exp = 0;
    mag_set_fmpz_2exp_fmpz_lower(z, x, exp);
}

MAG_INLINE void
mag_mul_ui(mag_t z, const mag_t x, ulong y)
{
    mag_t t;
    mag_init(t);
    mag_set_ui(t, y);
    mag_mul(z, x, t);
    mag_clear(t);
}

MAG_INLINE void
mag_mul_ui_lower(mag_t z, const mag_t x, ulong y)
{
    mag_t t;
    mag_init(t);
    mag_set_ui_lower(t, y);
    mag_mul_lower(z, x, t);
    mag_clear(t);
}

MAG_INLINE void
mag_mul_fmpz(mag_t z, const mag_t x, const fmpz_t y)
{
    mag_t t;
    mag_init(t);
    mag_set_fmpz(t, y);
    mag_mul(z, x, t);
    mag_clear(t);
}

MAG_INLINE void
mag_mul_fmpz_lower(mag_t z, const mag_t x, const fmpz_t y)
{
    mag_t t;
    mag_init(t);
    mag_set_fmpz_lower(t, y);
    mag_mul_lower(z, x, t);
    mag_clear(t);
}

MAG_INLINE void
mag_div_ui(mag_t z, const mag_t x, ulong y)
{
    mag_t t;
    mag_init(t);
    mag_set_ui_lower(t, y);
    mag_div(z, x, t);
    mag_clear(t);
}

MAG_INLINE void
mag_div_fmpz(mag_t z, const mag_t x, const fmpz_t y)
{
    mag_t t;
    mag_init(t);
    mag_set_fmpz_lower(t, y);
    mag_div(z, x, t);
    mag_clear(t);
}

#ifdef __cplusplus
}
#endif

#endif

