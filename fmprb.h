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
        flint_printf("ABORT: nan or negative radius!\n");
        abort();
    }
#endif

    fmpr_clear(fmprb_midref(x));
    fmpr_clear(fmprb_radref(x));
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

void fmprb_set_round(fmprb_t z, const fmprb_t x, slong prec);

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
fmprb_neg_round(fmprb_t x, const fmprb_t y, slong prec)
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
fmprb_set_si(fmprb_t x, slong y)
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
fmprb_set_round_fmpz_2exp(fmprb_t y, const fmpz_t x, const fmpz_t exp, slong prec)
{
    slong r = fmpr_set_round_fmpz_2exp(fmprb_midref(y), x, exp, prec, FMPR_RND_DOWN);
    fmpr_set_error_result(fmprb_radref(y), fmprb_midref(y), r);
}

static __inline__ void
fmprb_set_round_fmpz(fmprb_t y, const fmpz_t x, slong prec)
{
    slong r = fmpr_set_round_fmpz(fmprb_midref(y), x, prec, FMPR_RND_DOWN);
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
    fmprb_set_ui(f, UWORD(1));
}

static __inline__ void
fmprb_adjust(fmprb_t x)
{
    if (!fmprb_is_exact(x))
    {
        /* reduce precision here */
    }
}

void fmprb_add(fmprb_t z, const fmprb_t x, const fmprb_t y, slong prec);
void fmprb_add_ui(fmprb_t z, const fmprb_t x, ulong y, slong prec);
void fmprb_add_si(fmprb_t z, const fmprb_t x, slong y, slong prec);
void fmprb_add_fmpz(fmprb_t z, const fmprb_t x, const fmpz_t y, slong prec);
void fmprb_add_fmpr(fmprb_t z, const fmprb_t x, const fmpr_t y, slong prec);

void fmprb_addmul(fmprb_t z, const fmprb_t x, const fmprb_t y, slong prec);
void fmprb_addmul_ui(fmprb_t z, const fmprb_t x, ulong y, slong prec);
void fmprb_addmul_si(fmprb_t z, const fmprb_t x, slong y, slong prec);
void fmprb_addmul_fmpz(fmprb_t z, const fmprb_t x, const fmpz_t y, slong prec);

void fmprb_div(fmprb_t z, const fmprb_t x, const fmprb_t y, slong prec);
void fmprb_div_ui(fmprb_t z, const fmprb_t x, ulong y, slong prec);
void fmprb_div_si(fmprb_t z, const fmprb_t x, slong y, slong prec);
void fmprb_div_fmpz(fmprb_t z, const fmprb_t x, const fmpz_t y, slong prec);
void fmprb_div_fmpr(fmprb_t z, const fmprb_t x, const fmpr_t y, slong prec);
void fmprb_fmpz_div_fmpz(fmprb_t y, const fmpz_t num, const fmpz_t den, slong prec);
void fmprb_ui_div(fmprb_t z, ulong x, const fmprb_t y, slong prec);

static __inline__ void
fmprb_inv(fmprb_t y, const fmprb_t x, slong prec)
{
    fmprb_ui_div(y, 1, x, prec);
}

void fmprb_mul_fmpr_naive(fmprb_t z, const fmprb_t x, const fmpr_t y, slong prec);
void fmprb_mul_main_naive(fmprb_t z, const fmprb_t x, const fmprb_t y, slong prec);
void fmprb_mul_naive(fmprb_t z, const fmprb_t x, const fmprb_t y, slong prec);

void fmprb_mul(fmprb_t z, const fmprb_t x, const fmprb_t y, slong prec);
void fmprb_mul_ui(fmprb_t z, const fmprb_t x, ulong y, slong prec);
void fmprb_mul_si(fmprb_t z, const fmprb_t x, slong y, slong prec);
void fmprb_mul_fmpz(fmprb_t z, const fmprb_t x, const fmpz_t y, slong prec);
void fmprb_mul_fmpr(fmprb_t z, const fmprb_t x, const fmpr_t y, slong prec);

void fmprb_sqrt(fmprb_t z, const fmprb_t x, slong prec);
void fmprb_sqrt_ui(fmprb_t z, ulong x, slong prec);
void fmprb_sqrt_fmpz(fmprb_t z, const fmpz_t x, slong prec);

void fmprb_rsqrt(fmprb_t z, const fmprb_t x, slong prec);
void fmprb_rsqrt_ui(fmprb_t z, ulong x, slong prec);

void fmprb_sqrtpos(fmprb_t z, const fmprb_t x, slong prec);

void fmprb_hypot(fmprb_t z, const fmprb_t x, const fmprb_t y, slong prec);


void fmprb_sub(fmprb_t z, const fmprb_t x, const fmprb_t y, slong prec);
void fmprb_sub_ui(fmprb_t z, const fmprb_t x, ulong y, slong prec);
void fmprb_sub_si(fmprb_t z, const fmprb_t x, slong y, slong prec);
void fmprb_sub_fmpz(fmprb_t z, const fmprb_t x, const fmpz_t y, slong prec);

void fmprb_submul(fmprb_t z, const fmprb_t x, const fmprb_t y, slong prec);
void fmprb_submul_ui(fmprb_t z, const fmprb_t x, ulong y, slong prec);
void fmprb_submul_si(fmprb_t z, const fmprb_t x, slong y, slong prec);
void fmprb_submul_fmpz(fmprb_t z, const fmprb_t x, const fmpz_t y, slong prec);

static __inline__ void
fmprb_print(const fmprb_t x)
{
    fmpr_print(fmprb_midref(x));
    flint_printf(" +/- ");
    fmpr_print(fmprb_radref(x));
}

static __inline__ void
fmprb_printd(const fmprb_t x, slong digits)
{
    fmpr_printd(fmprb_midref(x), FLINT_ABS(digits));
    if (digits > 0)
    {
        flint_printf(" +/- ");
        fmpr_printd(fmprb_radref(x), 5);
    }
}

static __inline__ void
fmprb_mul_2exp_si(fmprb_t y, const fmprb_t x, slong e)
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
fmprb_set_fmpq(fmprb_t y, const fmpq_t x, slong prec)
{
    fmprb_fmpz_div_fmpz(y, fmpq_numref(x), fmpq_denref(x), prec);
}

int fmprb_contains_fmpr(const fmprb_t x, const fmpr_t y);
int fmprb_contains_fmpq(const fmprb_t x, const fmpq_t y);
int fmprb_contains_fmpz(const fmprb_t x, const fmpz_t y);
int fmprb_contains_si(const fmprb_t x, slong y);
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
fmprb_get_abs_ubound_fmpr(fmpr_t u, const fmprb_t x, slong prec)
{
    if (fmpr_sgn(fmprb_midref(x)) < 0)
        fmpr_sub(u, fmprb_midref(x), fmprb_radref(x), prec, FMPR_RND_UP);
    else
        fmpr_add(u, fmprb_midref(x), fmprb_radref(x), prec, FMPR_RND_UP);

    fmpr_abs(u, u);
}

static __inline__ void
fmprb_get_abs_lbound_fmpr(fmpr_t u, const fmprb_t x, slong prec)
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

void fmprb_set_interval_fmpr(fmprb_t x, const fmpr_t a, const fmpr_t b, slong prec);

static __inline__ slong
fmprb_bits(const fmprb_t x)
{
    return fmpr_bits(fmprb_midref(x));
}

void fmprb_add_error_fmpr(fmprb_t x, const fmpr_t err);
void fmprb_add_error_2exp_si(fmprb_t x, slong err);
void fmprb_add_error_2exp_fmpz(fmprb_t x, const fmpz_t err);
void fmprb_add_error(fmprb_t x, const fmprb_t error);

void fmprb_randtest(fmprb_t x, flint_rand_t state, slong prec, slong mag_bits);
void fmprb_randtest_exact(fmprb_t x, flint_rand_t state, slong prec, slong mag_bits);
void fmprb_randtest_wide(fmprb_t x, flint_rand_t state, slong prec, slong mag_bits);
void fmprb_randtest_precise(fmprb_t x, flint_rand_t state, slong prec, slong mag_bits);
void fmprb_randtest_special(fmprb_t x, flint_rand_t state, slong prec, slong mag_bits);

void fmprb_get_rand_fmpq(fmpq_t q, flint_rand_t state, const fmprb_t x, slong bits);

#ifdef __cplusplus
}
#endif

#endif
