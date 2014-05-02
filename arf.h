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

    2x2 mul code taken from MPFR 2.3.0
    (Copyright (C) 1991-2007 Free Software Foundation, Inc.)

******************************************************************************/

#ifndef ARF_H
#define ARF_H

#include <math.h>
#include "flint.h"
#include "fmpr.h"

#ifdef __cplusplus
extern "C" {
#endif

#define LIMB_ONE ((mp_limb_t) 1)
#define LIMB_ONES (-(mp_limb_t) 1)
#define LIMB_TOP (((mp_limb_t) 1) << (FLINT_BITS - 1))

#define MASK_LIMB(n, c) ((n) & (LIMB_ONES << (c)))

#define arf_rnd_t fmpr_rnd_t
#define ARF_RND_DOWN FMPR_RND_DOWN
#define ARF_RND_UP FMPR_RND_UP
#define ARF_RND_FLOOR FMPR_RND_FLOOR
#define ARF_RND_CEIL FMPR_RND_CEIL
#define ARF_RND_NEAR FMPR_RND_NEAR

static __inline__ int
arf_rounds_down(arf_rnd_t rnd, int sgnbit)
{
    if (rnd == ARF_RND_DOWN) return 1;
    if (rnd == ARF_RND_UP) return 0;
    if (rnd == ARF_RND_FLOOR) return !sgnbit;
    return sgnbit;
}

static __inline__ int
arf_rounds_up(arf_rnd_t rnd, int sgnbit)
{
    if (rnd == ARF_RND_DOWN) return 0;
    if (rnd == ARF_RND_UP) return 1;
    if (rnd == ARF_RND_FLOOR) return sgnbit;
    return !sgnbit;
}

static __inline__ mpfr_rnd_t
arf_rnd_to_mpfr(arf_rnd_t rnd)
{
    if (rnd == ARF_RND_DOWN) return MPFR_RNDZ;
    else if (rnd == ARF_RND_UP) return MPFR_RNDA;
    else if (rnd == ARF_RND_FLOOR) return MPFR_RNDD;
    else if (rnd == ARF_RND_CEIL) return MPFR_RNDU;
    else return MPFR_RNDN;
}

/* Allow 'infinite' precision for operations where a result can be computed exactly. */
#define ARF_PREC_EXACT LONG_MAX

#define ARF_PREC_ADD(prec,extra) ((prec) == ARF_PREC_EXACT ? ARF_PREC_EXACT : (prec) + (extra))

#define ARF_RESULT_EXACT 0
#define ARF_RESULT_INEXACT 1

/* Range where we can skip fmpz overflow checks for exponent manipulation. */
#define ARF_MAX_LAGOM_EXP (COEFF_MAX / 4)
#define ARF_MIN_LAGOM_EXP (-ARF_MAX_LAGOM_EXP)

/* Exponent values used to encode special values. */
#define ARF_EXP_ZERO 0
#define ARF_EXP_NAN COEFF_MIN
#define ARF_EXP_POS_INF (COEFF_MIN+1)
#define ARF_EXP_NEG_INF (COEFF_MIN+2)

/* Direct access to the exponent. */
#define ARF_EXP(x) ((x)->exp)
#define ARF_EXPREF(x) (&(x)->exp)

/* Finite and with lagom big exponents. */
#define ARF_IS_LAGOM(x) (ARF_EXP(x) >= ARF_MIN_LAGOM_EXP && \
                         ARF_EXP(x) <= ARF_MAX_LAGOM_EXP)

/* More than two limbs (needs pointer). */
#define ARF_HAS_PTR(x) ((x)->size > ((2 << 1) + 1))

/* Raw size field (encodes both limb size and sign). */
#define ARF_XSIZE(x) ((x)->size)

/* Construct size field value from size in limbs and sign bit. */
#define ARF_MAKE_XSIZE(size, sgnbit) ((((mp_size_t) size) << 1) | sgnbit)

/* The limb size, and the sign bit. */
#define ARF_SIZE(x) (ARF_XSIZE(x) >> 1)
#define ARF_SGNBIT(x) (ARF_XSIZE(x) & 1)

/* Assumes non-special value */
#define ARF_NEG(x) (ARF_XSIZE(x) ^= 1)

/* Note: may also be hardcoded in a few places. */
#define ARF_NOPTR_LIMBS 2

/* Direct access to the limb data. */
#define ARF_NOPTR_D(x)   ((x)->d.noptr.d)
#define ARF_PTR_D(x)     ((x)->d.ptr.d)
#define ARF_PTR_ALLOC(x) ((x)->d.ptr.alloc)

/* Encoding for special values. */
#define ARF_IS_SPECIAL(x) (ARF_XSIZE(x) == 0)

/* Value is +/- a power of two */
#define ARF_IS_POW2(x) (ARF_SIZE(x) == 1) && (ARF_NOPTR_D(x)[0] == LIMB_TOP)

/* To set a special value, first call this and then set the exponent. */
#define ARF_MAKE_SPECIAL(x)         \
    do {                            \
        fmpz_clear(ARF_EXPREF(x));  \
        ARF_DEMOTE(x);              \
        ARF_XSIZE(x) = 0;           \
    } while (0)


typedef struct
{
    mp_limb_t d[ARF_NOPTR_LIMBS];
}
mantissa_noptr_struct;

typedef struct
{
    mp_size_t alloc;
    mp_ptr d;
}
mantissa_ptr_struct;

typedef union
{
    mantissa_noptr_struct noptr;
    mantissa_ptr_struct ptr;
}
mantissa_struct;

typedef struct
{
    fmpz exp;
    mp_size_t size;
    mantissa_struct d;
}
arf_struct;

typedef arf_struct arf_t[1];
typedef arf_struct * arf_ptr;
typedef const arf_struct * arf_srcptr;

/* Warning: does not set size! -- also doesn't demote exponent. */
#define ARF_DEMOTE(x)                 \
    do {                              \
        if (ARF_HAS_PTR(x))           \
            flint_free(ARF_PTR_D(x)); \
    } while (0)

/* Get mpn pointer and size (xptr, xn) for read-only use. */
#define ARF_GET_MPN_READONLY(xptr, xn, x)   \
    do {                                    \
        xn = ARF_SIZE(x);                   \
        if (xn <= ARF_NOPTR_LIMBS)          \
            xptr = ARF_NOPTR_D(x);          \
        else                                \
            xptr = ARF_PTR_D(x);            \
    } while (0)

/* Get mpn pointer xptr for writing *exactly* xn limbs to x. */
#define ARF_GET_MPN_WRITE(xptr, xn, x)                      \
    do {                                                    \
        mp_size_t __xn = (xn);                              \
        if ((__xn) <= ARF_NOPTR_LIMBS)                      \
        {                                                   \
            ARF_DEMOTE(x);                                  \
            xptr = ARF_NOPTR_D(x);                          \
        }                                                   \
        else                                                \
        {                                                   \
            if (!ARF_HAS_PTR(x))                            \
            {                                               \
                ARF_PTR_D(x) = flint_malloc((__xn) *        \
                        sizeof(mp_limb_t));                 \
                ARF_PTR_ALLOC(x) = (__xn);                  \
            }                                               \
            else if (ARF_PTR_ALLOC(x) < (__xn))             \
            {                                               \
                ARF_PTR_D(x) = flint_realloc(ARF_PTR_D(x),  \
                        (xn) * sizeof(mp_limb_t));          \
                ARF_PTR_ALLOC(x) = (__xn);                  \
            }                                               \
            xptr = ARF_PTR_D(x);                            \
        }                                                   \
        ARF_XSIZE(x) = ARF_MAKE_XSIZE(__xn, 0);             \
    } while (0)

static __inline__ void
arf_init(arf_t x)
{
    fmpz_init(ARF_EXPREF(x));
    ARF_XSIZE(x) = 0;
}

static __inline__ void
arf_clear(arf_t x)
{
    fmpz_clear(ARF_EXPREF(x));
    ARF_DEMOTE(x);
}

static __inline__ void
arf_zero(arf_t x)
{
    ARF_MAKE_SPECIAL(x);
    ARF_EXP(x) = ARF_EXP_ZERO;
}

static __inline__ void
arf_pos_inf(arf_t x)
{
    ARF_MAKE_SPECIAL(x);
    ARF_EXP(x) = ARF_EXP_POS_INF;
}

static __inline__ void
arf_neg_inf(arf_t x)
{
    ARF_MAKE_SPECIAL(x);
    ARF_EXP(x) = ARF_EXP_NEG_INF;
}

static __inline__ void
arf_nan(arf_t x)
{
    ARF_MAKE_SPECIAL(x);
    ARF_EXP(x) = ARF_EXP_NAN;
}

static __inline__ int
arf_is_special(const arf_t x)
{
    return ARF_IS_SPECIAL(x);
}

static __inline__ int
arf_is_zero(const arf_t x)
{
    return ARF_IS_SPECIAL(x) && (ARF_EXP(x) == ARF_EXP_ZERO);
}

static __inline__ int
arf_is_pos_inf(const arf_t x)
{
    return ARF_IS_SPECIAL(x) && (ARF_EXP(x) == ARF_EXP_POS_INF);
}

static __inline__ int
arf_is_neg_inf(const arf_t x)
{
    return ARF_IS_SPECIAL(x) && (ARF_EXP(x) == ARF_EXP_NEG_INF);
}

static __inline__ int
arf_is_nan(const arf_t x)
{
    return ARF_IS_SPECIAL(x) && (ARF_EXP(x) == ARF_EXP_NAN);
}

static __inline__ int
arf_is_normal(const arf_t x)
{
    return !ARF_IS_SPECIAL(x);
}

static __inline__ int
arf_is_finite(const arf_t x)
{
    return !ARF_IS_SPECIAL(x) || (ARF_EXP(x) == ARF_EXP_ZERO);
}

static __inline__ int
arf_is_inf(const arf_t x)
{
    return ARF_IS_SPECIAL(x) && (ARF_EXP(x) == ARF_EXP_POS_INF ||
                                 ARF_EXP(x) == ARF_EXP_NEG_INF);
}

static __inline__ void
arf_one(arf_t x)
{
    fmpz_clear(ARF_EXPREF(x));
    ARF_DEMOTE(x);
    ARF_EXP(x) = 1;
    ARF_XSIZE(x) = ARF_MAKE_XSIZE(1, 0);
    ARF_NOPTR_D(x)[0] = LIMB_TOP;
}

static __inline__ int
arf_is_one(const arf_t x)
{
    return (ARF_EXP(x) == 1) && (ARF_XSIZE(x) == ARF_MAKE_XSIZE(1, 0))
                             && ARF_NOPTR_D(x)[0] == LIMB_TOP;
}

static __inline__ int
arf_sgn(const arf_t x)
{
    if (arf_is_special(x))
    {
        if (arf_is_zero(x) || arf_is_nan(x))
            return 0;
        return arf_is_pos_inf(x) ? 1 : -1;
    }
    else
    {
        return ARF_SGNBIT(x) ? -1 : 1;
    }
}

int arf_cmp(const arf_t x, const arf_t y);

int arf_cmpabs(const arf_t x, const arf_t y);

static __inline__ void
arf_swap(arf_t y, arf_t x)
{
    if (x != y)
    {
        arf_struct t = *x;
        *x = *y;
        *y = t;
    }
}

static __inline__ void
arf_set(arf_t y, const arf_t x)
{
    if (x != y)
    {
        /* Fast path */
        if (!COEFF_IS_MPZ(ARF_EXP(x)) && !COEFF_IS_MPZ(ARF_EXP(y)))
            ARF_EXP(y) = ARF_EXP(x);
        else
            fmpz_set(ARF_EXPREF(y), ARF_EXPREF(x));

        /* Fast path */
        if (!ARF_HAS_PTR(x))
        {
            ARF_DEMOTE(y);
            (y)->d = (x)->d;
        }
        else
        {
            mp_ptr yptr;
            mp_srcptr xptr;
            mp_size_t n;

            ARF_GET_MPN_READONLY(xptr, n, x);
            ARF_GET_MPN_WRITE(yptr, n, y);
            flint_mpn_copyi(yptr, xptr, n);
        }

        /* Important. */
        ARF_XSIZE(y) = ARF_XSIZE(x);
    }
}

static __inline__ void
arf_neg(arf_t y, const arf_t x)
{
    arf_set(y, x);

    if (arf_is_special(y))
    {
        if (arf_is_pos_inf(y))
            arf_neg_inf(y);
        else if (arf_is_neg_inf(y))
            arf_pos_inf(y);
    }
    else
    {
        ARF_NEG(y);
    }
}

static __inline__ void
arf_set_ui(arf_t x, ulong v)
{
    ARF_DEMOTE(x);
    _fmpz_demote(ARF_EXPREF(x));

    if (v == 0)
    {
        ARF_EXP(x) = ARF_EXP_ZERO;
        ARF_XSIZE(x) = 0;
    }
    else
    {
        unsigned int c;
        count_leading_zeros(c, v);
        ARF_EXP(x) = FLINT_BITS - c;
        ARF_NOPTR_D(x)[0] = v << c;
        ARF_XSIZE(x) = ARF_MAKE_XSIZE(1, 0);
    }
}

static __inline__ void
arf_set_si(arf_t x, long v)
{
    arf_set_ui(x, FLINT_ABS(v));
    if (v < 0)
        ARF_NEG(x);
}

/* Assumes xn > 0, x[0] != 0. */
/* TBD: 1, 2 limb versions */
void arf_set_mpn(arf_t y, mp_srcptr x, mp_size_t xn, int sgnbit);

static __inline__ void
arf_set_mpz(arf_t y, const mpz_t x)
{
    long size = x->_mp_size;

    if (size == 0)
        arf_zero(y);
    else
        arf_set_mpn(y, x->_mp_d, FLINT_ABS(size), size < 0);
}

static __inline__ void
arf_set_fmpz(arf_t y, const fmpz_t x)
{
    if (!COEFF_IS_MPZ(*x))
        arf_set_si(y, *x);
    else
        arf_set_mpz(y, COEFF_TO_PTR(*x));
}

int _arf_set_round_ui(arf_t x, ulong v, int sgnbit, long prec, arf_rnd_t rnd);

/* Assumes xn > 0, x[0] != 0. */
/* TBD: 1, 2 limb versions, top-aligned version */
int
_arf_set_round_mpn(arf_t y, long * exp_shift, mp_srcptr x, mp_size_t xn,
    int sgnbit, long prec, arf_rnd_t rnd);

static __inline__ int
arf_set_round_ui(arf_t x, ulong v, long prec, arf_rnd_t rnd)
{
    return _arf_set_round_ui(x, v, 0, prec, rnd);
}

static __inline__ int
arf_set_round_si(arf_t x, long v, long prec, arf_rnd_t rnd)
{
    return _arf_set_round_ui(x, FLINT_ABS(v), v < 0, prec, rnd);
}

static __inline__ int
arf_set_round_mpz(arf_t y, const mpz_t x, long prec, arf_rnd_t rnd)
{
    int inexact;
    long size = x->_mp_size;
    long fix;

    if (size == 0)
    {
        arf_zero(y);
        return 0;
    }

    inexact = _arf_set_round_mpn(y, &fix, x->_mp_d, FLINT_ABS(size),
        (size < 0), prec, rnd);
    _fmpz_demote(ARF_EXPREF(y));
    ARF_EXP(y) = FLINT_ABS(size) * FLINT_BITS + fix;
    return inexact;
}

static __inline__ int
arf_set_round_fmpz(arf_t y, const fmpz_t x, long prec, arf_rnd_t rnd)
{
    if (!COEFF_IS_MPZ(*x))
        return arf_set_round_si(y, *x, prec, rnd);
    else
        return arf_set_round_mpz(y, COEFF_TO_PTR(*x), prec, rnd);
}

int arf_set_round(arf_t y, const arf_t x, long prec, arf_rnd_t rnd);

int arf_neg_round(arf_t y, const arf_t x, long prec, arf_rnd_t rnd);

void arf_get_fmpr(fmpr_t y, const arf_t x);

void arf_set_fmpr(arf_t y, const fmpr_t x);

int arf_get_mpfr(mpfr_t x, const arf_t y, mpfr_rnd_t rnd);

void arf_set_mpfr(arf_t x, const mpfr_t y);

int arf_equal(const arf_t x, const arf_t y);

static __inline__ void
arf_min(arf_t z, const arf_t a, const arf_t b)
{
    if (arf_cmp(a, b) <= 0)
        arf_set(z, a);
    else
        arf_set(z, b);
}

static __inline__ void
arf_max(arf_t z, const arf_t a, const arf_t b)
{
    if (arf_cmp(a, b) > 0)
        arf_set(z, a);
    else
        arf_set(z, b);
}

static __inline__ void
arf_abs(arf_t y, const arf_t x)
{
    if (arf_sgn(x) < 0)
        arf_neg(y, x);
    else
        arf_set(y, x);
}

static __inline__ long
arf_bits(const arf_t x)
{
    if (arf_is_special(x))
        return 0;
    else
    {
        mp_srcptr xp;
        mp_size_t xn;
        long c;

        ARF_GET_MPN_READONLY(xp, xn, x);
        count_trailing_zeros(c, xp[0]);
        return xn * FLINT_BITS - c;
    }
}

static __inline__ void
arf_bot(fmpz_t e, const arf_t x)
{
    if (arf_is_special(x))
        fmpz_zero(e);
    else
        fmpz_sub_si(e, ARF_EXPREF(x), arf_bits(x));
}

static __inline__ int
arf_is_int(const arf_t x)
{
    if (arf_is_special(x))
        return arf_is_zero(x);
    else
    {
        fmpz_t t;
        int r;
        fmpz_init(t);
        arf_bot(t, x);
        r = fmpz_sgn(t) >= 0;
        fmpz_clear(t);
        return r;
    }
}

static __inline__ int
arf_is_int_2exp_si(const arf_t x, long e)
{
    if (arf_is_special(x))
        return arf_is_zero(x);
    else
    {
        fmpz_t t;
        int r;
        fmpz_init(t);
        arf_bot(t, x);
        r = fmpz_cmp_si(t, e) >= 0;
        fmpz_clear(t);
        return r;
    }
}

int arf_cmp_2exp_si(const arf_t x, long e);

int arf_cmpabs_2exp_si(const arf_t x, long e);

static __inline__ void
arf_set_si_2exp_si(arf_t x, long man, long exp)
{
    arf_set_si(x, man);
    if (man != 0)
        fmpz_add_si_inline(ARF_EXPREF(x), ARF_EXPREF(x), exp);
}

static __inline__ void
arf_set_ui_2exp_si(arf_t x, ulong man, long exp)
{
    arf_set_ui(x, man);
    if (man != 0)
        fmpz_add_si_inline(ARF_EXPREF(x), ARF_EXPREF(x), exp);
}

static __inline__ void
arf_mul_2exp_si(arf_t y, const arf_t x, long e)
{
    arf_set(y, x);
    if (!arf_is_special(y))
        fmpz_add_si_inline(ARF_EXPREF(y), ARF_EXPREF(y), e);
}

static __inline__ void
arf_mul_2exp_fmpz(arf_t y, const arf_t x, const fmpz_t e)
{
    arf_set(y, x);
    if (!arf_is_special(y))
        fmpz_add_inline(ARF_EXPREF(y), ARF_EXPREF(y), e);
}

static __inline__ int
arf_set_round_fmpz_2exp(arf_t y, const fmpz_t x, const fmpz_t exp, long prec, arf_rnd_t rnd)
{
    if (fmpz_is_zero(x))
    {
        arf_zero(y);
        return 0;
    }
    else
    {
        int r = arf_set_round_fmpz(y, x, prec, rnd);
        fmpz_add_inline(ARF_EXPREF(y), ARF_EXPREF(y), exp);
        return r;
    }
}

static __inline__ void
arf_abs_bound_lt_2exp_fmpz(fmpz_t b, const arf_t x)
{
    if (arf_is_special(x))
        fmpz_zero(b);
    else
        fmpz_set(b, ARF_EXPREF(x));
}

static __inline__ void
arf_abs_bound_le_2exp_fmpz(fmpz_t b, const arf_t x)
{
    if (arf_is_special(x))
        fmpz_zero(b);
    else if (ARF_IS_POW2(x))
        fmpz_sub_ui(b, ARF_EXPREF(x), 1);
    else
        fmpz_set(b, ARF_EXPREF(x));
}

long arf_abs_bound_lt_2exp_si(const arf_t x);

void arf_get_fmpz_2exp(fmpz_t man, fmpz_t exp, const arf_t x);

static __inline__ void
arf_set_fmpz_2exp(arf_t x, const fmpz_t man, const fmpz_t exp)
{
    arf_set_fmpz(x, man);
    if (!arf_is_zero(x))
        fmpz_add_inline(ARF_EXPREF(x), ARF_EXPREF(x), exp);
}

void arf_debug(const arf_t x);

void arf_print(const arf_t x);

void arf_randtest(arf_t x, flint_rand_t state, long bits, long mag_bits);

void arf_randtest_not_zero(arf_t x, flint_rand_t state, long bits, long mag_bits);

void arf_randtest_special(arf_t x, flint_rand_t state, long bits, long mag_bits);

#define ADD2_FAST_MAX (COEFF_MAX / 4)
#define ADD2_FAST_MIN (-ADD2_FAST_MAX)

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

#define MUL_MPFR_MIN_LIMBS 25
#define MUL_MPFR_MAX_LIMBS 10000

#define nn_mul_2x1(r2, r1, r0, a1, a0, b0)                  \
    do {                                                    \
        mp_limb_t t1;                                       \
        umul_ppmm(r1, r0, a0, b0);                          \
        umul_ppmm(r2, t1, a1, b0);                          \
        add_ssaaaa(r2, r1, r2, r1, 0, t1);                  \
    } while (0)

#define nn_mul_2x2(r3, r2, r1, r0, a1, a0, b1, b0)          \
    do {                                                    \
        mp_limb_t t1, t2, t3;                               \
        umul_ppmm(r1, r0, a0, b0);                          \
        umul_ppmm(r2, t1, a1, b0);                          \
        add_ssaaaa(r2, r1, r2, r1, 0, t1);                  \
        umul_ppmm(t1, t2, a0, b1);                          \
        umul_ppmm(r3, t3, a1, b1);                          \
        add_ssaaaa(r3, t1, r3, t1, 0, t3);                  \
        add_ssaaaa(r2, r1, r2, r1, t1, t2);                 \
        r3 += r2 < t1;                                      \
    } while (0)

#define ARF_MPN_MUL(_z, _x, _xn, _y, _yn) \
    if ((_xn) == (_yn)) \
    { \
        if ((_xn) == 1) \
        { \
            umul_ppmm((_z)[1], (_z)[0], (_x)[0], (_y)[0]); \
        } \
        else if ((_xn) == 2) \
        { \
            mp_limb_t __x1, __x0, __y1, __y0; \
            __x0 = (_x)[0]; \
            __x1 = (_x)[1]; \
            __y0 = (_y)[0]; \
            __y1 = (_y)[1]; \
            nn_mul_2x2((_z)[3], (_z)[2], (_z)[1], (_z)[0], __x1, __x0, __y1, __y0); \
        } \
        else if ((_x) == (_y)) \
        { \
            mpn_sqr((_z), (_x), (_xn)); \
        } \
        else \
        { \
            mpn_mul_n((_z), (_x), (_y), (_xn)); \
        } \
    } \
    else if ((_xn) > (_yn)) \
    { \
        if ((_yn) == 1) \
            (_z)[(_xn) + (_yn) - 1] = mpn_mul_1((_z), (_x), (_xn), (_y)[0]); \
        else \
            mpn_mul((_z), (_x), (_xn), (_y), (_yn)); \
    } \
    else \
    { \
        if ((_xn) == 1) \
            (_z)[(_xn) + (_yn) - 1] = mpn_mul_1((_z), (_y), (_yn), (_x)[0]); \
        else \
            mpn_mul((_z), (_y), (_yn), (_x), (_xn)); \
    }

#define ARF_MUL_STACK_ALLOC 40
#define ARF_MUL_TLS_ALLOC 1000

extern TLS_PREFIX mp_ptr __arf_mul_tmp;
extern TLS_PREFIX long __arf_mul_alloc;

extern void _arf_mul_tmp_cleanup(void);

#define ARF_MUL_TMP_DECL \
    mp_limb_t tmp_stack[ARF_MUL_STACK_ALLOC]; \

#define ARF_MUL_TMP_ALLOC(tmp, alloc) \
    if (alloc <= ARF_MUL_STACK_ALLOC) \
    { \
        tmp = tmp_stack; \
    } \
    else if (alloc <= ARF_MUL_TLS_ALLOC) \
    { \
        if (__arf_mul_alloc < alloc) \
        { \
            if (__arf_mul_alloc == 0) \
            { \
                flint_register_cleanup_function(_arf_mul_tmp_cleanup); \
            } \
            __arf_mul_tmp = flint_realloc(__arf_mul_tmp, sizeof(mp_limb_t) * alloc); \
            __arf_mul_alloc = alloc; \
        } \
        tmp = __arf_mul_tmp; \
    } \
    else \
    { \
        tmp = flint_malloc(sizeof(mp_limb_t) * alloc); \
    }

#define ARF_MUL_TMP_FREE(tmp, alloc) \
    if (alloc > ARF_MUL_TLS_ALLOC) \
        flint_free(tmp);

void arf_mul_special(arf_t z, const arf_t x, const arf_t y);

int arf_mul_via_mpfr(arf_t z, const arf_t x, const arf_t y, long prec, arf_rnd_t rnd);

int arf_mul_rnd_any(arf_ptr z, arf_srcptr x, arf_srcptr y, long prec, arf_rnd_t rnd);

int arf_mul_rnd_down(arf_ptr z, arf_srcptr x, arf_srcptr y, long prec);

#define arf_mul(z, x, y, prec, rnd)              \
    ((rnd == FMPR_RND_DOWN)                      \
        ? arf_mul_rnd_down(z, x, y, prec)        \
        : arf_mul_rnd_any(z, x, y, prec, rnd))

#define ARF_ADD_STACK_ALLOC 40
#define ARF_ADD_TLS_ALLOC 1000

extern TLS_PREFIX mp_ptr __arf_add_tmp;
extern TLS_PREFIX long __arf_add_alloc;

extern void _arf_add_tmp_cleanup(void);

#define ARF_ADD_TMP_DECL \
    mp_limb_t tmp_stack[ARF_ADD_STACK_ALLOC]; \

#define ARF_ADD_TMP_ALLOC(tmp, alloc) \
    if (alloc <= ARF_ADD_STACK_ALLOC) \
    { \
        tmp = tmp_stack; \
    } \
    else if (alloc <= ARF_ADD_TLS_ALLOC) \
    { \
        if (__arf_add_alloc < alloc) \
        { \
            if (__arf_add_alloc == 0) \
            { \
                flint_register_cleanup_function(_arf_add_tmp_cleanup); \
            } \
            __arf_add_tmp = flint_realloc(__arf_add_tmp, sizeof(mp_limb_t) * alloc); \
            __arf_add_alloc = alloc; \
        } \
        tmp = __arf_add_tmp; \
    } \
    else \
    { \
        tmp = flint_malloc(sizeof(mp_limb_t) * alloc); \
    }

#define ARF_ADD_TMP_FREE(tmp, alloc) \
    if (alloc > ARF_ADD_TLS_ALLOC) \
        flint_free(tmp);

int _arf_add_mpn(arf_t z, mp_srcptr xp, mp_size_t xn, int xsgnbit,
    const fmpz_t xexp, mp_srcptr yp, mp_size_t yn, int ysgnbit,
    mp_bitcnt_t shift, long prec, arf_rnd_t rnd);

int arf_add(arf_ptr z, arf_srcptr x, arf_srcptr y, long prec, arf_rnd_t rnd);
int arf_add_si(arf_ptr z, arf_srcptr x, long y, long prec, arf_rnd_t rnd);
int arf_add_ui(arf_ptr z, arf_srcptr x, ulong y, long prec, arf_rnd_t rnd);
int arf_add_fmpz(arf_ptr z, arf_srcptr x, const fmpz_t y, long prec, arf_rnd_t rnd);

int arf_sub(arf_ptr z, arf_srcptr x, arf_srcptr y, long prec, arf_rnd_t rnd);

#ifdef __cplusplus
}
#endif

#endif

