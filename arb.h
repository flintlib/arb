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

#ifndef ARB_H
#define ARB_H

#include "fmprb.h"
#include "mag.h"
#include "arf.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    arf_struct mid;
    mag_struct rad;
}
arb_struct;

typedef arb_struct arb_t[1];
typedef arb_struct * arb_ptr;
typedef const arb_struct * arb_srcptr;

#define ARB_MIDREF(x) (&(x)->mid)
#define ARB_RADREF(x) (&(x)->rad)

#define arb_midref ARB_MIDREF
#define arb_radref ARB_RADREF

#define ARB_IS_LAGOM(x) (ARF_IS_LAGOM(ARB_MIDREF(x)) && MAG_IS_LAGOM(ARB_RADREF(x)))

#define ARB_RND ARF_RND_DOWN

static __inline__ void
arb_init(arb_t x)
{
    arf_init(ARB_MIDREF(x));
    mag_init(ARB_RADREF(x));
}

static __inline__ void
arb_clear(arb_t x)
{
    arf_clear(ARB_MIDREF(x));
    mag_clear(ARB_RADREF(x));
}

static __inline__ arb_ptr
_arb_vec_init(long n)
{
    long i;
    arb_ptr v = (arb_ptr) flint_malloc(sizeof(arb_struct) * n);

    for (i = 0; i < n; i++)
        arb_init(v + i);

    return v;
}

static __inline__ void
_arb_vec_clear(arb_ptr v, long n)
{
    long i;
    for (i = 0; i < n; i++)
        arb_clear(v + i);
    flint_free(v);
}

static __inline__ void
arb_set_fmprb(arb_t x, const fmprb_t y)
{
    arf_set_fmpr(ARB_MIDREF(x), fmprb_midref(y));
    mag_set_fmpr(ARB_RADREF(x), fmprb_radref(y));
}

static __inline__ void
arb_get_fmprb(fmprb_t x, const arb_t y)
{
    arf_get_fmpr(fmprb_midref(x), ARB_MIDREF(y));
    mag_get_fmpr(fmprb_radref(x), ARB_RADREF(y));
}

static __inline__ int
arb_is_exact(const arb_t x)
{
    return mag_is_zero(arb_radref(x));
}

static __inline__ int
arb_equal(const arb_t x, const arb_t y)
{
    return arf_equal(arb_midref(x), arb_midref(y)) &&
           mag_equal(arb_radref(x), arb_radref(y));
}

static __inline__ void
arb_zero(arb_t x)
{
    arf_zero(arb_midref(x));
    mag_zero(arb_radref(x));
}

static __inline__ int
arb_is_zero(const arb_t x)
{
    return arf_is_zero(arb_midref(x)) && mag_is_zero(arb_radref(x));
}

static __inline__ void
arb_pos_inf(arb_t x)
{
    arf_pos_inf(arb_midref(x));
    mag_zero(arb_radref(x));
}

static __inline__ void
arb_neg_inf(arb_t x)
{
    arf_neg_inf(arb_midref(x));
    mag_zero(arb_radref(x));
}

static __inline__ void
arb_zero_pm_inf(arb_t x)
{
    arf_zero(arb_midref(x));
    mag_inf(arb_radref(x));
}

static __inline__ void
arb_indeterminate(arb_t x)
{
    arf_nan(arb_midref(x));
    mag_inf(arb_radref(x));
}

static __inline__ int
arb_is_finite(const arb_t x)
{
    return arf_is_finite(arb_midref(x)) && mag_is_finite(arb_radref(x));
}

static __inline__ void
arb_set(arb_t x, const arb_t y)
{
    arf_set(arb_midref(x), arb_midref(y));
    mag_set(arb_radref(x), arb_radref(y));
}

static __inline__ void
arb_swap(arb_t x, arb_t y)
{
    arb_struct t = *x;
    *x = *y;
    *y = t;
}

void arb_set_round(arb_t z, const arb_t x, long prec);

void arb_trim(arb_t y, const arb_t x);

static __inline__ void
arb_neg(arb_t x, const arb_t y)
{
    arf_neg(arb_midref(x), arb_midref(y));
    mag_set(arb_radref(x), arb_radref(y));
}

static __inline__ void
arb_neg_round(arb_t x, const arb_t y, long prec)
{
    arb_set_round(x, y, prec);
    arb_neg(x, x);
}

static __inline__ void
arb_abs(arb_t x, const arb_t y)
{
    arf_abs(arb_midref(x), arb_midref(y));
    mag_set(arb_radref(x), arb_radref(y));
}

static __inline__ void
arb_set_arf(arb_t x, const arf_t y)
{
    arf_set(arb_midref(x), y);
    mag_zero(arb_radref(x));
}

static __inline__ void
arb_set_si(arb_t x, long y)
{
    arf_set_si(arb_midref(x), y);
    mag_zero(arb_radref(x));
}

static __inline__ void
arb_set_ui(arb_t x, ulong y)
{
    arf_set_ui(arb_midref(x), y);
    mag_zero(arb_radref(x));
}

static __inline__ void
arb_set_fmpz(arb_t x, const fmpz_t y)
{
    arf_set_fmpz(arb_midref(x), y);
    mag_zero(arb_radref(x));
}

static __inline__ void
arb_set_fmpz_2exp(arb_t x, const fmpz_t y, const fmpz_t exp)
{
    arf_set_fmpz_2exp(arb_midref(x), y, exp);
    mag_zero(arb_radref(x));
}

static __inline__ void
arb_set_round_fmpz_2exp(arb_t y, const fmpz_t x, const fmpz_t exp, long prec)
{
    int inexact;
    inexact = arf_set_round_fmpz_2exp(arb_midref(y), x, exp, prec, ARB_RND);

    if (inexact)
        arf_mag_set_ulp(arb_radref(y), arb_midref(y), prec);
    else
        mag_zero(arb_radref(y));
}

static __inline__ void
arb_set_round_fmpz(arb_t y, const fmpz_t x, long prec)
{
    int inexact;
    inexact = arf_set_round_fmpz(arb_midref(y), x, prec, ARB_RND);

    if (inexact)
        arf_mag_set_ulp(arb_radref(y), arb_midref(y), prec);
    else
        mag_zero(arb_radref(y));
}

static __inline__ int
arb_is_one(const arb_t f)
{
    return arf_is_one(arb_midref(f)) && mag_is_zero(arb_radref(f));
}

static __inline__ void
arb_one(arb_t f)
{
    arb_set_ui(f, 1UL);
}

static __inline__ void
arb_print(const arb_t x)
{
    arf_print(arb_midref(x));
    printf(" +/- ");
    mag_print(arb_radref(x));
}

static __inline__ void
arb_printd(const arb_t x, long digits)
{
    arf_printd(arb_midref(x), FLINT_ABS(digits));
    if (digits > 0)
    {
        printf(" +/- ");
        mag_printd(arb_radref(x), 5);
    }
}

static __inline__ void
arb_mul_2exp_si(arb_t y, const arb_t x, long e)
{
    arf_mul_2exp_si(arb_midref(y), arb_midref(x), e);
    mag_mul_2exp_si(arb_radref(y), arb_radref(x), e);
}

static __inline__ void
arb_mul_2exp_fmpz(arb_t y, const arb_t x, const fmpz_t e)
{
    arf_mul_2exp_fmpz(arb_midref(y), arb_midref(x), e);
    mag_mul_2exp_fmpz(arb_radref(y), arb_radref(x), e);
}

static __inline__ int
arb_is_int(const arb_t x)
{
    return mag_is_zero(arb_radref(x)) &&
           arf_is_int(arb_midref(x));
}

static __inline__ int
arb_contains_zero(const arb_t x)
{
    return arf_cmpabs_mag(arb_midref(x), arb_radref(x)) <= 0;
}

static __inline__ int
arb_is_nonzero(const arb_t x)
{
    return !arb_contains_zero(x);
}

static __inline__ int
arb_is_positive(const arb_t x)
{
    return (arf_sgn(arb_midref(x)) > 0) &&
        (arf_mag_cmpabs(arb_radref(x), arb_midref(x)) < 0) &&
         !arf_is_nan(arb_midref(x));
}

static __inline__ int
arb_is_nonnegative(const arb_t x)
{
    return (arf_sgn(arb_midref(x)) >= 0) &&
        (arf_mag_cmpabs(arb_radref(x), arb_midref(x)) <= 0) &&
         !arf_is_nan(arb_midref(x));
}

static __inline__ int
arb_is_negative(const arb_t x)
{
    return (arf_sgn(arb_midref(x)) < 0) &&
        (arf_mag_cmpabs(arb_radref(x), arb_midref(x)) < 0) &&
         !arf_is_nan(arb_midref(x));
}

static __inline__ int
arb_is_nonpositive(const arb_t x)
{
    return (arf_sgn(arb_midref(x)) <= 0) &&
        (arf_mag_cmpabs(arb_radref(x), arb_midref(x)) <= 0) &&
         !arf_is_nan(arb_midref(x));
}

static __inline__ int
arb_contains_negative(const arb_t x)
{
    return (arf_sgn(arb_midref(x)) < 0) ||
        (arf_mag_cmpabs(arb_radref(x), arb_midref(x)) > 0)
        || arf_is_nan(arb_midref(x));
}

static __inline__ int
arb_contains_nonpositive(const arb_t x)
{
    return (arf_sgn(arb_midref(x)) <= 0) ||
        (arf_mag_cmpabs(arb_radref(x), arb_midref(x)) >= 0)
        || arf_is_nan(arb_midref(x));
}

static __inline__ int
arb_contains_positive(const arb_t x)
{
    return (arf_sgn(arb_midref(x)) > 0) ||
        (arf_mag_cmpabs(arb_radref(x), arb_midref(x)) > 0)
        || arf_is_nan(arb_midref(x));
}

static __inline__ int
arb_contains_nonnegative(const arb_t x)
{
    return (arf_sgn(arb_midref(x)) >= 0) ||
        (arf_mag_cmpabs(arb_radref(x), arb_midref(x)) >= 0)
        || arf_is_nan(arb_midref(x));
}

/* TODO: get/set mag versions */

static __inline__ void
arb_get_abs_ubound_arf(arf_t u, const arb_t x, long prec)
{
    arf_t t;
    arf_init_set_mag_shallow(t, arb_radref(x));

    if (arf_sgn(arb_midref(x)) < 0)
        arf_sub(u, arb_midref(x), t, prec, ARF_RND_UP);
    else
        arf_add(u, arb_midref(x), t, prec, ARF_RND_UP);

    arf_abs(u, u);
}

static __inline__ void
arb_get_abs_lbound_fmpr(arf_t u, const arb_t x, long prec)
{
    arf_t t;
    arf_init_set_mag_shallow(t, arb_radref(x));

    if (arf_sgn(arb_midref(x)) > 0)
    {
        arf_sub(u, arb_midref(x), t, prec, ARF_RND_UP);
    }
    else
    {
        arf_add(u, arb_midref(x), t, prec, ARF_RND_UP);
        arf_neg(u, u);
    }

    if (arf_sgn(u) < 0)
        arf_zero(u);
}

static __inline__ long
arb_rel_error_bits(const arb_t x)
{
    fmpz_t midmag, radmag;
    arf_t t;
    long result;

    if (mag_is_zero(arb_radref(x)))
        return -ARF_PREC_EXACT;
    if (arf_is_special(arb_midref(x)) || mag_is_inf(arb_radref(x)))
        return ARF_PREC_EXACT;

    fmpz_init(midmag);
    fmpz_init(radmag);

    arf_abs_bound_lt_2exp_fmpz(midmag, arb_midref(x));

    arf_init_set_mag_shallow(t, arb_radref(x)); /* no need to free */
    arf_abs_bound_lt_2exp_fmpz(radmag, t);
    fmpz_add_ui(radmag, radmag, 1);

    result = _fmpz_sub_small(radmag, midmag);

    fmpz_clear(midmag);
    fmpz_clear(radmag);

    return result;
}

static __inline__ long
arb_rel_accuracy_bits(const arb_t x)
{
    return -arb_rel_error_bits(x);
}

static __inline__ long
arb_bits(const arb_t x)
{
    return arf_bits(arb_midref(x));
}


/* vector functions */

static __inline__ void
_arb_vec_zero(arb_ptr A, long n)
{
    long i;
    for (i = 0; i < n; i++)
        arb_zero(A + i);
}

static __inline__ int
_arb_vec_is_zero(arb_srcptr vec, long len)
{
    long i;
    for (i = 0; i < len; i++)
        if (!arb_is_zero(vec + i))
            return 0;
    return 1;
}

static __inline__ void
_arb_vec_set(arb_ptr res, arb_srcptr vec, long len)
{
    long i;
    for (i = 0; i < len; i++)
        arb_set(res + i, vec + i);
}

static __inline__ void
_arb_vec_set_round(arb_ptr res, arb_srcptr vec, long len, long prec)
{
    long i;
    for (i = 0; i < len; i++)
        arb_set_round(res + i, vec + i, prec);
}

static __inline__ void
_arb_vec_swap(arb_ptr res, arb_ptr vec, long len)
{
    long i;
    for (i = 0; i < len; i++)
        arb_swap(res + i, vec + i);
}

static __inline__ void
_arb_vec_neg(arb_ptr B, arb_srcptr A, long n)
{
    long i;
    for (i = 0; i < n; i++)
        arb_neg(B + i, A + i);
}

/*
static __inline__ void
_arb_vec_sub(arb_ptr C, arb_srcptr A,
    arb_srcptr B, long n, long prec)
{
    long i;
    for (i = 0; i < n; i++)
        arb_sub(C + i, A + i, B + i, prec);
}

static __inline__ void
_arb_vec_add(arb_ptr C, arb_srcptr A,
    arb_srcptr B, long n, long prec)
{
    long i;
    for (i = 0; i < n; i++)
        arb_add(C + i, A + i, B + i, prec);
}

static __inline__ void
_arb_vec_scalar_mul(arb_ptr res, arb_srcptr vec,
    long len, const arb_t c, long prec)
{
    long i;
    for (i = 0; i < len; i++)
        arb_mul(res + i, vec + i, c, prec);
}

static __inline__ void
_arb_vec_scalar_div(arb_ptr res, arb_srcptr vec,
    long len, const arb_t c, long prec)
{
    long i;
    for (i = 0; i < len; i++)
        arb_div(res + i, vec + i, c, prec);
}

static __inline__ void
_arb_vec_scalar_mul_fmpz(arb_ptr res, arb_srcptr vec,
    long len, const fmpz_t c, long prec)
{
    long i;
    arf_t t;
    arf_init(t);
    arf_set_fmpz(t, c);
    for (i = 0; i < len; i++)
        arb_mul_arf(res + i, vec + i, t, prec);
    arf_clear(t);
}

static __inline__ void
_arb_vec_scalar_mul_2exp_si(arb_ptr res, arb_srcptr src, long len, long c)
{
    long i;
    for (i = 0; i < len; i++)
        arb_mul_2exp_si(res + i, src + i, c);
}

static __inline__ void
_arb_vec_scalar_addmul(arb_ptr res, arb_srcptr vec,
    long len, const arb_t c, long prec)
{
    long i;
    for (i = 0; i < len; i++)
        arb_addmul(res + i, vec + i, c, prec);
}
*/

/* TODO: mag version? */
static __inline__ void
_arb_vec_get_abs_ubound_arf(arf_t bound, arb_srcptr vec,
        long len, long prec)
{
    arf_t t;
    long i;

    if (len < 1)
    {
        arf_zero(bound);
    }
    else
    {
        arb_get_abs_ubound_arf(bound, vec, prec);
        arf_init(t);
        for (i = 1; i < len; i++)
        {
            arb_get_abs_ubound_arf(t, vec + i, prec);
            if (arf_cmp(t, bound) > 0)
                arf_set(bound, t);
        }
        arf_clear(t);
    }
}

static __inline__ long
_arb_vec_bits(arb_srcptr x, long len)
{
    long i, b, c;

    b = 0;
    for (i = 0; i < len; i++)
    {
        c = arb_bits(x + i);
        b = FLINT_MAX(b, c);
    }

    return b;
}

/*
static __inline__ void
_arb_vec_set_powers(arb_ptr xs, const arb_t x, long len, long prec)
{
    long i;

    for (i = 0; i < len; i++)
    {
        if (i == 0)
            arb_one(xs + i);
        else if (i == 1)
            arb_set_round(xs + i, x, prec);
        else if (i % 2 == 0)
            arb_mul(xs + i, xs + i / 2, xs + i / 2, prec);
        else
            arb_mul(xs + i, xs + i - 1, x, prec);
    }
}
*/

/* TODO: mag verision ? */
/*
static __inline__ void
_arb_vec_add_error_arf_vec(arb_ptr res, arf_srcptr err, long len)
{
    long i;
    for (i = 0; i < len; i++)
        arb_add_error_arf(res + i, err + i);
}
*/

static __inline__ void
_arb_vec_indeterminate(arb_ptr vec, long len)
{
    long i;
    for (i = 0; i < len; i++)
        arb_indeterminate(vec + i);
}

static __inline__ void
_arb_vec_trim(arb_ptr res, arb_srcptr vec, long len)
{
    long i;
    for (i = 0; i < len; i++)
        arb_trim(res + i, vec + i);
}

void arb_mul(arb_t z, const arb_t x, const arb_t y, long prec);

void arb_addmul(arb_t z, const arb_t x, const arb_t y, long prec);

/*
static __inline__ void
arb_set_fmpq(arb_t y, const fmpq_t x, long prec)
{
    arb_fmpz_div_fmpz(y, fmpq_numref(x), fmpq_denref(x), prec);
}
*/

void arb_randtest_exact(arb_t x, flint_rand_t state, long prec, long mag_bits);

void arb_randtest_wide(arb_t x, flint_rand_t state, long prec, long mag_bits);

void arb_randtest_precise(arb_t x, flint_rand_t state, long prec, long mag_bits);

void arb_randtest(arb_t x, flint_rand_t state, long prec, long mag_bits);

void arb_randtest_special(arb_t x, flint_rand_t state, long prec, long mag_bits);

void arb_add_error_arf(arb_t x, const arf_t err);

void arb_add_error_2exp_si(arb_t x, long err);

void arb_add_error_2exp_fmpz(arb_t x, const fmpz_t err);

void arb_add_error(arb_t x, const arb_t error);

int arb_contains_arf(const arb_t x, const arf_t y);

int arb_contains_fmpq(const arb_t x, const fmpq_t y);

int arb_contains_fmpz(const arb_t x, const fmpz_t y);

int arb_contains_si(const arb_t x, long y);

int arb_contains_mpfr(const arb_t x, const mpfr_t y);

int arb_overlaps(const arb_t x, const arb_t y);

int arb_contains(const arb_t x, const arb_t y);

void arb_get_interval_fmpz_2exp(fmpz_t a, fmpz_t b, fmpz_t exp, const arb_t x);
int arb_get_unique_fmpz(fmpz_t z, const arb_t x);
void arb_set_interval_arf(arb_t x, const arf_t a, const arf_t b, long prec);
void arb_union(arb_t z, const arb_t x, const arb_t y, long prec);
void arb_get_rand_fmpq(fmpq_t q, flint_rand_t state, const arb_t x, long bits);

void arb_add(arb_t z, const arb_t x, const arb_t y, long prec);
void arb_add_arf(arb_t z, const arb_t x, const arf_t y, long prec);
void arb_add_ui(arb_t z, const arb_t x, ulong y, long prec);
void arb_add_si(arb_t z, const arb_t x, long y, long prec);
void arb_add_fmpz(arb_t z, const arb_t x, const fmpz_t y, long prec);

void arb_sub(arb_t z, const arb_t x, const arb_t y, long prec);
void arb_sub_arf(arb_t z, const arb_t x, const arf_t y, long prec);
void arb_sub_ui(arb_t z, const arb_t x, ulong y, long prec);
void arb_sub_si(arb_t z, const arb_t x, long y, long prec);
void arb_sub_fmpz(arb_t z, const arb_t x, const fmpz_t y, long prec);

#ifdef __cplusplus
}
#endif

#endif

