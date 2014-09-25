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

#ifndef FMPRB_H
#define FMPRB_H

#include "fmpr.h"
#include "fmpz_poly.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    fmpr_struct mid;
    fmpr_struct rad;
}
fmprb_struct;

typedef fmprb_struct fmprb_t[1];
typedef fmprb_struct * fmprb_ptr;
typedef const fmprb_struct * fmprb_srcptr;

#define fmprb_midref(x) (&(x)->mid)
#define fmprb_radref(x) (&(x)->rad)


#define FMPRB_RAD_PREC 30
#define FMPRB_MIN_ADJUSTED_PREC 

static __inline__ void
fmprb_init(fmprb_t x)
{
    fmpr_init(fmprb_midref(x));
    fmpr_init(fmprb_radref(x));
}

#define FMPRB_DEBUG 0

static __inline__ void
fmprb_clear(fmprb_t x)
{
#if FMPRB_DEBUG
    if (fmpr_is_nan(fmprb_radref(x)) || fmpr_sgn(fmprb_radref(x)) < 0)
    {
        printf("ABORT: nan or negative radius!\n");
        abort();
    }
#endif

    fmpr_clear(fmprb_midref(x));
    fmpr_clear(fmprb_radref(x));
}

static __inline__ fmprb_ptr
_fmprb_vec_init(long n)
{
    long i;
    fmprb_ptr v = (fmprb_ptr) flint_malloc(sizeof(fmprb_struct) * n);

    for (i = 0; i < n; i++)
        fmprb_init(v + i);

    return v;
}

static __inline__ void
_fmprb_vec_clear(fmprb_ptr v, long n)
{
    long i;
    for (i = 0; i < n; i++)
        fmprb_clear(v + i);
    flint_free(v);
}

static __inline__ int
fmprb_is_exact(const fmprb_t x)
{
    return fmpr_is_zero(fmprb_radref(x));
}

static __inline__ int
fmprb_equal(const fmprb_t x, const fmprb_t y)
{
    return fmpr_equal(fmprb_midref(x), fmprb_midref(y)) &&
            fmpr_equal(fmprb_radref(x), fmprb_radref(y));
}

static __inline__ void
fmprb_zero(fmprb_t x)
{
    fmpr_zero(fmprb_midref(x));
    fmpr_zero(fmprb_radref(x));
}

static __inline__ int
fmprb_is_zero(const fmprb_t x)
{
    return fmpr_is_zero(fmprb_midref(x)) && fmpr_is_zero(fmprb_radref(x));
}

static __inline__ void
fmprb_pos_inf(fmprb_t x)
{
    fmpr_pos_inf(fmprb_midref(x));
    fmpr_zero(fmprb_radref(x));
}

static __inline__ void
fmprb_neg_inf(fmprb_t x)
{
    fmpr_neg_inf(fmprb_midref(x));
    fmpr_zero(fmprb_radref(x));
}

static __inline__ void
fmprb_zero_pm_inf(fmprb_t x)
{
    fmpr_zero(fmprb_midref(x));
    fmpr_pos_inf(fmprb_radref(x));
}

static __inline__ void
fmprb_indeterminate(fmprb_t x)
{
    fmpr_nan(fmprb_midref(x));
    fmpr_pos_inf(fmprb_radref(x));
}

static __inline__ int
fmprb_is_finite(const fmprb_t x)
{
    return fmpr_is_finite(fmprb_midref(x)) && fmpr_is_finite(fmprb_radref(x));
}

static __inline__ void
fmprb_set(fmprb_t x, const fmprb_t y)
{
    fmpr_set(fmprb_midref(x), fmprb_midref(y));
    fmpr_set(fmprb_radref(x), fmprb_radref(y));
}

void fmprb_set_round(fmprb_t z, const fmprb_t x, long prec);

void fmprb_trim(fmprb_t y, const fmprb_t x);

static __inline__ void
fmprb_swap(fmprb_t x, fmprb_t y)
{
    fmprb_struct t = *x;
    *x = *y;
    *y = t;
}

static __inline__ void
fmprb_neg(fmprb_t x, const fmprb_t y)
{
    fmpr_neg(fmprb_midref(x), fmprb_midref(y));
    fmpr_set(fmprb_radref(x), fmprb_radref(y));
}

static __inline__ void
fmprb_neg_round(fmprb_t x, const fmprb_t y, long prec)
{
    fmprb_set_round(x, y, prec);
    fmprb_neg(x, x);
}

static __inline__ void
fmprb_abs(fmprb_t x, const fmprb_t y)
{
    fmpr_abs(fmprb_midref(x), fmprb_midref(y));
    fmpr_set(fmprb_radref(x), fmprb_radref(y));
}

static __inline__ void
fmprb_set_fmpr(fmprb_t x, const fmpr_t y)
{
    fmpr_set(fmprb_midref(x), y);
    fmpr_zero(fmprb_radref(x));
}

static __inline__ void
fmprb_set_si(fmprb_t x, long y)
{
    fmpr_set_si(fmprb_midref(x), y);
    fmpr_zero(fmprb_radref(x));
}

static __inline__ void
fmprb_set_ui(fmprb_t x, ulong y)
{
    fmpr_set_ui(fmprb_midref(x), y);
    fmpr_zero(fmprb_radref(x));
}

static __inline__ void
fmprb_set_fmpz(fmprb_t x, const fmpz_t y)
{
    fmpr_set_fmpz(fmprb_midref(x), y);
    fmpr_zero(fmprb_radref(x));
}

static __inline__ void
fmprb_set_fmpz_2exp(fmprb_t x, const fmpz_t y, const fmpz_t exp)
{
    fmpr_set_fmpz_2exp(fmprb_midref(x), y, exp);
    fmpr_zero(fmprb_radref(x));
}


static __inline__ void
fmprb_set_round_fmpz_2exp(fmprb_t y, const fmpz_t x, const fmpz_t exp, long prec)
{
    long r = fmpr_set_round_fmpz_2exp(fmprb_midref(y), x, exp, prec, FMPR_RND_DOWN);
    fmpr_set_error_result(fmprb_radref(y), fmprb_midref(y), r);
}

static __inline__ void
fmprb_set_round_fmpz(fmprb_t y, const fmpz_t x, long prec)
{
    long r = fmpr_set_round_fmpz(fmprb_midref(y), x, prec, FMPR_RND_DOWN);
    fmpr_set_error_result(fmprb_radref(y), fmprb_midref(y), r);
}

static __inline__ int
fmprb_is_one(const fmprb_t f)
{
    return fmpr_is_one(fmprb_midref(f)) && fmpr_is_zero(fmprb_radref(f));
}

static __inline__ void
fmprb_one(fmprb_t f)
{
    fmprb_set_ui(f, 1UL);
}

static __inline__ void
fmprb_adjust(fmprb_t x)
{
    if (!fmprb_is_exact(x))
    {
        /* reduce precision here */
    }
}

void fmprb_add(fmprb_t z, const fmprb_t x, const fmprb_t y, long prec);
void fmprb_add_ui(fmprb_t z, const fmprb_t x, ulong y, long prec);
void fmprb_add_si(fmprb_t z, const fmprb_t x, long y, long prec);
void fmprb_add_fmpz(fmprb_t z, const fmprb_t x, const fmpz_t y, long prec);
void fmprb_add_fmpr(fmprb_t z, const fmprb_t x, const fmpr_t y, long prec);

void fmprb_addmul(fmprb_t z, const fmprb_t x, const fmprb_t y, long prec);
void fmprb_addmul_ui(fmprb_t z, const fmprb_t x, ulong y, long prec);
void fmprb_addmul_si(fmprb_t z, const fmprb_t x, long y, long prec);
void fmprb_addmul_fmpz(fmprb_t z, const fmprb_t x, const fmpz_t y, long prec);

void fmprb_div(fmprb_t z, const fmprb_t x, const fmprb_t y, long prec);
void fmprb_div_ui(fmprb_t z, const fmprb_t x, ulong y, long prec);
void fmprb_div_si(fmprb_t z, const fmprb_t x, long y, long prec);
void fmprb_div_fmpz(fmprb_t z, const fmprb_t x, const fmpz_t y, long prec);
void fmprb_div_fmpr(fmprb_t z, const fmprb_t x, const fmpr_t y, long prec);
void fmprb_fmpz_div_fmpz(fmprb_t y, const fmpz_t num, const fmpz_t den, long prec);
void fmprb_ui_div(fmprb_t z, ulong x, const fmprb_t y, long prec);

static __inline__ void
fmprb_inv(fmprb_t y, const fmprb_t x, long prec)
{
    fmprb_ui_div(y, 1, x, prec);
}

void fmprb_div_2expm1_ui(fmprb_t y, const fmprb_t x, ulong n, long prec);


void fmprb_mul_fmpr_naive(fmprb_t z, const fmprb_t x, const fmpr_t y, long prec);
void fmprb_mul_main_naive(fmprb_t z, const fmprb_t x, const fmprb_t y, long prec);
void fmprb_mul_naive(fmprb_t z, const fmprb_t x, const fmprb_t y, long prec);

void fmprb_mul(fmprb_t z, const fmprb_t x, const fmprb_t y, long prec);
void fmprb_mul_ui(fmprb_t z, const fmprb_t x, ulong y, long prec);
void fmprb_mul_si(fmprb_t z, const fmprb_t x, long y, long prec);
void fmprb_mul_fmpz(fmprb_t z, const fmprb_t x, const fmpz_t y, long prec);
void fmprb_mul_fmpr(fmprb_t z, const fmprb_t x, const fmpr_t y, long prec);

void fmprb_sqrt(fmprb_t z, const fmprb_t x, long prec);
void fmprb_sqrt_ui(fmprb_t z, ulong x, long prec);
void fmprb_sqrt_fmpz(fmprb_t z, const fmpz_t x, long prec);

void fmprb_rsqrt(fmprb_t z, const fmprb_t x, long prec);
void fmprb_rsqrt_ui(fmprb_t z, ulong x, long prec);

void fmprb_sqrtpos(fmprb_t z, const fmprb_t x, long prec);

void fmprb_hypot(fmprb_t z, const fmprb_t x, const fmprb_t y, long prec);


void fmprb_sub(fmprb_t z, const fmprb_t x, const fmprb_t y, long prec);
void fmprb_sub_ui(fmprb_t z, const fmprb_t x, ulong y, long prec);
void fmprb_sub_si(fmprb_t z, const fmprb_t x, long y, long prec);
void fmprb_sub_fmpz(fmprb_t z, const fmprb_t x, const fmpz_t y, long prec);

void fmprb_submul(fmprb_t z, const fmprb_t x, const fmprb_t y, long prec);
void fmprb_submul_ui(fmprb_t z, const fmprb_t x, ulong y, long prec);
void fmprb_submul_si(fmprb_t z, const fmprb_t x, long y, long prec);
void fmprb_submul_fmpz(fmprb_t z, const fmprb_t x, const fmpz_t y, long prec);

void fmprb_root(fmprb_t z, const fmprb_t x, ulong k, long prec);

void fmprb_agm(fmprb_t z, const fmprb_t x, const fmprb_t y, long prec);

static __inline__ void
fmprb_print(const fmprb_t x)
{
    fmpr_print(fmprb_midref(x));
    printf(" +/- ");
    fmpr_print(fmprb_radref(x));
}

static __inline__ void
fmprb_printd(const fmprb_t x, long digits)
{
    fmpr_printd(fmprb_midref(x), FLINT_ABS(digits));
    if (digits > 0)
    {
        printf(" +/- ");
        fmpr_printd(fmprb_radref(x), 5);
    }
}

static __inline__ void
fmprb_mul_2exp_si(fmprb_t y, const fmprb_t x, long e)
{
    fmpr_mul_2exp_si(fmprb_midref(y), fmprb_midref(x), e);
    fmpr_mul_2exp_si(fmprb_radref(y), fmprb_radref(x), e);
}

static __inline__ void
fmprb_mul_2exp_fmpz(fmprb_t y, const fmprb_t x, const fmpz_t e)
{
    fmpr_mul_2exp_fmpz(fmprb_midref(y), fmprb_midref(x), e);
    fmpr_mul_2exp_fmpz(fmprb_radref(y), fmprb_radref(x), e);
}

static __inline__ void
fmprb_set_fmpq(fmprb_t y, const fmpq_t x, long prec)
{
    fmprb_fmpz_div_fmpz(y, fmpq_numref(x), fmpq_denref(x), prec);
}

int fmprb_contains_fmpr(const fmprb_t x, const fmpr_t y);
int fmprb_contains_fmpq(const fmprb_t x, const fmpq_t y);
int fmprb_contains_fmpz(const fmprb_t x, const fmpz_t y);
int fmprb_contains_si(const fmprb_t x, long y);
int fmprb_contains_mpfr(const fmprb_t x, const mpfr_t y);
int fmprb_contains_zero(const fmprb_t x);

int fmprb_overlaps(const fmprb_t x, const fmprb_t y);

int fmprb_contains(const fmprb_t x, const fmprb_t y);

static __inline__ int
fmprb_is_int(const fmprb_t x)
{
    return fmpr_is_zero(fmprb_radref(x)) &&
           fmpr_is_int(fmprb_midref(x));
}

static __inline__ int
fmprb_is_nonzero(const fmprb_t x)
{
    return !fmprb_contains_zero(x);
}

static __inline__ int
fmprb_is_positive(const fmprb_t x)
{
    return (fmpr_sgn(fmprb_midref(x)) > 0) &&
        (fmpr_cmpabs(fmprb_radref(x), fmprb_midref(x)) < 0) &&
         !fmpr_is_nan(fmprb_midref(x));
}

static __inline__ int
fmprb_is_nonnegative(const fmprb_t x)
{
    return (fmpr_sgn(fmprb_midref(x)) >= 0) &&
        (fmpr_cmpabs(fmprb_radref(x), fmprb_midref(x)) <= 0) &&
         !fmpr_is_nan(fmprb_midref(x));
}

static __inline__ int
fmprb_is_negative(const fmprb_t x)
{
    return (fmpr_sgn(fmprb_midref(x)) < 0) &&
        (fmpr_cmpabs(fmprb_radref(x), fmprb_midref(x)) < 0) &&
         !fmpr_is_nan(fmprb_midref(x));
}

static __inline__ int
fmprb_is_nonpositive(const fmprb_t x)
{
    return (fmpr_sgn(fmprb_midref(x)) <= 0) &&
        (fmpr_cmpabs(fmprb_radref(x), fmprb_midref(x)) <= 0) &&
         !fmpr_is_nan(fmprb_midref(x));
}

static __inline__ int
fmprb_contains_negative(const fmprb_t x)
{
    return (fmpr_sgn(fmprb_midref(x)) < 0) ||
        (fmpr_cmpabs(fmprb_radref(x), fmprb_midref(x)) > 0)
        || fmpr_is_nan(fmprb_midref(x));
}

static __inline__ int
fmprb_contains_nonpositive(const fmprb_t x)
{
    return (fmpr_sgn(fmprb_midref(x)) <= 0) ||
        (fmpr_cmpabs(fmprb_radref(x), fmprb_midref(x)) >= 0)
        || fmpr_is_nan(fmprb_midref(x));
}

static __inline__ int
fmprb_contains_positive(const fmprb_t x)
{
    return (fmpr_sgn(fmprb_midref(x)) > 0) ||
        (fmpr_cmpabs(fmprb_radref(x), fmprb_midref(x)) > 0)
        || fmpr_is_nan(fmprb_midref(x));
}

static __inline__ int
fmprb_contains_nonnegative(const fmprb_t x)
{
    return (fmpr_sgn(fmprb_midref(x)) >= 0) ||
        (fmpr_cmpabs(fmprb_radref(x), fmprb_midref(x)) >= 0)
        || fmpr_is_nan(fmprb_midref(x));
}

static __inline__ void
fmprb_get_abs_ubound_fmpr(fmpr_t u, const fmprb_t x, long prec)
{
    if (fmpr_sgn(fmprb_midref(x)) < 0)
        fmpr_sub(u, fmprb_midref(x), fmprb_radref(x), prec, FMPR_RND_UP);
    else
        fmpr_add(u, fmprb_midref(x), fmprb_radref(x), prec, FMPR_RND_UP);

    fmpr_abs(u, u);
}

static __inline__ void
fmprb_get_abs_lbound_fmpr(fmpr_t u, const fmprb_t x, long prec)
{
    if (fmpr_sgn(fmprb_midref(x)) > 0)
    {
        fmpr_sub(u, fmprb_midref(x), fmprb_radref(x), prec, FMPR_RND_DOWN);
    }
    else
    {
        fmpr_add(u, fmprb_midref(x), fmprb_radref(x), prec, FMPR_RND_DOWN);
        fmpr_neg(u, u);
    }

    if (fmpr_sgn(u) < 0)
        fmpr_zero(u);
}

void fmprb_get_interval_fmpz_2exp(fmpz_t a, fmpz_t b, fmpz_t exp, const fmprb_t x);

int fmprb_get_unique_fmpz(fmpz_t z, const fmprb_t x);

void fmprb_set_interval_fmpr(fmprb_t x, const fmpr_t a, const fmpr_t b, long prec);

void fmprb_union(fmprb_t z, const fmprb_t x, const fmprb_t y, long prec);

static __inline__ long
fmprb_rel_error_bits(const fmprb_t x)
{
    fmpz_t midmag, radmag;
    long result;

    if (fmpr_is_zero(fmprb_radref(x)))
        return -FMPR_PREC_EXACT;
    if (fmpr_is_special(fmprb_midref(x)) || fmpr_is_special(fmprb_radref(x)))
        return FMPR_PREC_EXACT;

    fmpz_init(midmag);
    fmpz_init(radmag);

    fmpr_abs_bound_lt_2exp_fmpz(midmag, fmprb_midref(x));
    fmpr_abs_bound_lt_2exp_fmpz(radmag, fmprb_radref(x));
    fmpz_add_ui(radmag, radmag, 1);

    result = _fmpz_sub_small(radmag, midmag);

    fmpz_clear(midmag);
    fmpz_clear(radmag);

    return result;
}

static __inline__ long
fmprb_rel_accuracy_bits(const fmprb_t x)
{
    return -fmprb_rel_error_bits(x);
}

static __inline__ long
fmprb_bits(const fmprb_t x)
{
    return fmpr_bits(fmprb_midref(x));
}

void fmprb_add_error_fmpr(fmprb_t x, const fmpr_t err);
void fmprb_add_error_2exp_si(fmprb_t x, long err);
void fmprb_add_error_2exp_fmpz(fmprb_t x, const fmpz_t err);
void fmprb_add_error(fmprb_t x, const fmprb_t error);

void fmprb_randtest(fmprb_t x, flint_rand_t state, long prec, long mag_bits);
void fmprb_randtest_exact(fmprb_t x, flint_rand_t state, long prec, long mag_bits);
void fmprb_randtest_wide(fmprb_t x, flint_rand_t state, long prec, long mag_bits);
void fmprb_randtest_precise(fmprb_t x, flint_rand_t state, long prec, long mag_bits);
void fmprb_randtest_special(fmprb_t x, flint_rand_t state, long prec, long mag_bits);

void fmprb_get_rand_fmpq(fmpq_t q, flint_rand_t state, const fmprb_t x, long bits);

#define DEF_CACHED_CONSTANT(name, comp_func) \
    TLS_PREFIX long name ## _cached_prec = 0; \
    TLS_PREFIX fmprb_t name ## _cached_value; \
    void name ## _cleanup(void) \
    { \
        fmprb_clear(name ## _cached_value); \
        name ## _cached_prec = 0; \
    } \
    void name(fmprb_t x, long prec) \
    { \
        if (name ## _cached_prec < prec) \
        { \
            if (name ## _cached_prec == 0) \
            { \
                fmprb_init(name ## _cached_value); \
                flint_register_cleanup_function(name ## _cleanup); \
            } \
            comp_func(name ## _cached_value, prec); \
            name ## _cached_prec = prec; \
        } \
        fmprb_set_round(x, name ## _cached_value, prec); \
    }

/* vector functions */

static __inline__ void
_fmprb_vec_zero(fmprb_ptr A, long n)
{
    long i;
    for (i = 0; i < n; i++)
        fmprb_zero(A + i);
}

static __inline__ int
_fmprb_vec_is_zero(fmprb_srcptr vec, long len)
{
    long i;
    for (i = 0; i < len; i++)
        if (!fmprb_is_zero(vec + i))
            return 0;
    return 1;
}

static __inline__ void
_fmprb_vec_set(fmprb_ptr res, fmprb_srcptr vec, long len)
{
    long i;
    for (i = 0; i < len; i++)
        fmprb_set(res + i, vec + i);
}

static __inline__ void
_fmprb_vec_set_round(fmprb_ptr res, fmprb_srcptr vec, long len, long prec)
{
    long i;
    for (i = 0; i < len; i++)
        fmprb_set_round(res + i, vec + i, prec);
}

static __inline__ void
_fmprb_vec_swap(fmprb_ptr res, fmprb_ptr vec, long len)
{
    long i;
    for (i = 0; i < len; i++)
        fmprb_swap(res + i, vec + i);
}

static __inline__ void
_fmprb_vec_neg(fmprb_ptr B, fmprb_srcptr A, long n)
{
    long i;
    for (i = 0; i < n; i++)
        fmprb_neg(B + i, A + i);
}

static __inline__ void
_fmprb_vec_sub(fmprb_ptr C, fmprb_srcptr A,
    fmprb_srcptr B, long n, long prec)
{
    long i;
    for (i = 0; i < n; i++)
        fmprb_sub(C + i, A + i, B + i, prec);
}

static __inline__ void
_fmprb_vec_add(fmprb_ptr C, fmprb_srcptr A,
    fmprb_srcptr B, long n, long prec)
{
    long i;
    for (i = 0; i < n; i++)
        fmprb_add(C + i, A + i, B + i, prec);
}

static __inline__ void
_fmprb_vec_scalar_mul(fmprb_ptr res, fmprb_srcptr vec,
    long len, const fmprb_t c, long prec)
{
    long i;
    for (i = 0; i < len; i++)
        fmprb_mul(res + i, vec + i, c, prec);
}

static __inline__ void
_fmprb_vec_scalar_div(fmprb_ptr res, fmprb_srcptr vec,
    long len, const fmprb_t c, long prec)
{
    long i;
    for (i = 0; i < len; i++)
        fmprb_div(res + i, vec + i, c, prec);
}

static __inline__ void
_fmprb_vec_scalar_mul_fmpz(fmprb_ptr res, fmprb_srcptr vec,
    long len, const fmpz_t c, long prec)
{
    long i;
    fmpr_t t;
    fmpr_init(t);
    fmpr_set_fmpz(t, c);
    for (i = 0; i < len; i++)
        fmprb_mul_fmpr(res + i, vec + i, t, prec);
    fmpr_clear(t);
}

static __inline__ void
_fmprb_vec_scalar_mul_2exp_si(fmprb_ptr res, fmprb_srcptr src, long len, long c)
{
    long i;
    for (i = 0; i < len; i++)
        fmprb_mul_2exp_si(res + i, src + i, c);
}

static __inline__ void
_fmprb_vec_scalar_addmul(fmprb_ptr res, fmprb_srcptr vec,
    long len, const fmprb_t c, long prec)
{
    long i;
    for (i = 0; i < len; i++)
        fmprb_addmul(res + i, vec + i, c, prec);
}

static __inline__ void
_fmprb_vec_get_abs_ubound_fmpr(fmpr_t bound, fmprb_srcptr vec,
        long len, long prec)
{
    fmpr_t t;
    long i;

    if (len < 1)
    {
        fmpr_zero(bound);
    }
    else
    {
        fmprb_get_abs_ubound_fmpr(bound, vec, prec);
        fmpr_init(t);
        for (i = 1; i < len; i++)
        {
            fmprb_get_abs_ubound_fmpr(t, vec + i, prec);
            if (fmpr_cmp(t, bound) > 0)
                fmpr_set(bound, t);
        }
        fmpr_clear(t);
    }
}

static __inline__ long
_fmprb_vec_bits(fmprb_srcptr x, long len)
{
    long i, b, c;

    b = 0;
    for (i = 0; i < len; i++)
    {
        c = fmprb_bits(x + i);
        b = FLINT_MAX(b, c);
    }

    return b;
}

static __inline__ void
_fmprb_vec_set_powers(fmprb_ptr xs, const fmprb_t x, long len, long prec)
{
    long i;

    for (i = 0; i < len; i++)
    {
        if (i == 0)
            fmprb_one(xs + i);
        else if (i == 1)
            fmprb_set_round(xs + i, x, prec);
        else if (i % 2 == 0)
            fmprb_mul(xs + i, xs + i / 2, xs + i / 2, prec);
        else
            fmprb_mul(xs + i, xs + i - 1, x, prec);
    }
}

static __inline__ void
_fmprb_vec_add_error_fmpr_vec(fmprb_ptr res, fmpr_srcptr err, long len)
{
    long i;
    for (i = 0; i < len; i++)
        fmprb_add_error_fmpr(res + i, err + i);
}

static __inline__ void
_fmprb_vec_indeterminate(fmprb_ptr vec, long len)
{
    long i;
    for (i = 0; i < len; i++)
    {
        fmpr_nan(fmprb_midref(vec + i));
        fmpr_pos_inf(fmprb_radref(vec + i));
    }
}

static __inline__ void
_fmprb_vec_trim(fmprb_ptr res, fmprb_srcptr vec, long len)
{
    long i;
    for (i = 0; i < len; i++)
        fmprb_trim(res + i, vec + i);
}

#ifdef __cplusplus
}
#endif

#endif
