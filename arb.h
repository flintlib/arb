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

#ifdef ARB_INLINES_C
#define ARB_INLINE
#else
#define ARB_INLINE static __inline__
#endif

#include "fmprb.h"
#include "mag.h"
#include "arf.h"

#ifdef __cplusplus
extern "C" {
#endif

extern const char * arb_version;

typedef struct
{
    arf_struct mid;
    mag_struct rad;
}
arb_struct;

typedef arb_struct arb_t[1];
typedef arb_struct * arb_ptr;
typedef const arb_struct * arb_srcptr;

#define arb_midref(x) (&(x)->mid)
#define arb_radref(x) (&(x)->rad)

#define ARB_IS_LAGOM(x) (ARF_IS_LAGOM(arb_midref(x)) && MAG_IS_LAGOM(arb_radref(x)))

#define ARB_RND ARF_RND_DOWN

ARB_INLINE void
arb_init(arb_t x)
{
    arf_init(arb_midref(x));
    mag_init(arb_radref(x));
}

ARB_INLINE void
arb_clear(arb_t x)
{
    arf_clear(arb_midref(x));
    mag_clear(arb_radref(x));
}

ARB_INLINE arb_ptr
_arb_vec_init(long n)
{
    long i;
    arb_ptr v = (arb_ptr) flint_malloc(sizeof(arb_struct) * n);

    for (i = 0; i < n; i++)
        arb_init(v + i);

    return v;
}

ARB_INLINE void
_arb_vec_clear(arb_ptr v, long n)
{
    long i;
    for (i = 0; i < n; i++)
        arb_clear(v + i);
    flint_free(v);
}

ARB_INLINE void
arb_set_fmprb(arb_t x, const fmprb_t y)
{
    arf_set_fmpr(arb_midref(x), fmprb_midref(y));
    mag_set_fmpr(arb_radref(x), fmprb_radref(y));
}

ARB_INLINE void
arb_get_fmprb(fmprb_t x, const arb_t y)
{
    arf_get_fmpr(fmprb_midref(x), arb_midref(y));
    mag_get_fmpr(fmprb_radref(x), arb_radref(y));
}

ARB_INLINE int
arb_is_exact(const arb_t x)
{
    return mag_is_zero(arb_radref(x));
}

ARB_INLINE int
arb_equal(const arb_t x, const arb_t y)
{
    return arf_equal(arb_midref(x), arb_midref(y)) &&
           mag_equal(arb_radref(x), arb_radref(y));
}

ARB_INLINE void
arb_zero(arb_t x)
{
    arf_zero(arb_midref(x));
    mag_zero(arb_radref(x));
}

ARB_INLINE int
arb_is_zero(const arb_t x)
{
    return arf_is_zero(arb_midref(x)) && mag_is_zero(arb_radref(x));
}

ARB_INLINE void
arb_pos_inf(arb_t x)
{
    arf_pos_inf(arb_midref(x));
    mag_zero(arb_radref(x));
}

ARB_INLINE void
arb_neg_inf(arb_t x)
{
    arf_neg_inf(arb_midref(x));
    mag_zero(arb_radref(x));
}

ARB_INLINE void
arb_zero_pm_inf(arb_t x)
{
    arf_zero(arb_midref(x));
    mag_inf(arb_radref(x));
}

ARB_INLINE void
arb_indeterminate(arb_t x)
{
    arf_nan(arb_midref(x));
    mag_inf(arb_radref(x));
}

ARB_INLINE int
arb_is_finite(const arb_t x)
{
    return arf_is_finite(arb_midref(x)) && mag_is_finite(arb_radref(x));
}

void arb_set(arb_t x, const arb_t y);

ARB_INLINE void
arb_swap(arb_t x, arb_t y)
{
    arb_struct t = *x;
    *x = *y;
    *y = t;
}

void arb_set_round(arb_t z, const arb_t x, long prec);

void arb_trim(arb_t y, const arb_t x);

void arb_neg(arb_t x, const arb_t y);

void arb_neg_round(arb_t x, const arb_t y, long prec);

void arb_abs(arb_t x, const arb_t y);

ARB_INLINE void
arb_set_arf(arb_t x, const arf_t y)
{
    arf_set(arb_midref(x), y);
    mag_zero(arb_radref(x));
}

ARB_INLINE void
arb_set_si(arb_t x, long y)
{
    arf_set_si(arb_midref(x), y);
    mag_zero(arb_radref(x));
}

ARB_INLINE void
arb_set_ui(arb_t x, ulong y)
{
    arf_set_ui(arb_midref(x), y);
    mag_zero(arb_radref(x));
}

ARB_INLINE void
arb_set_fmpz(arb_t x, const fmpz_t y)
{
    arf_set_fmpz(arb_midref(x), y);
    mag_zero(arb_radref(x));
}

ARB_INLINE void
arb_set_fmpz_2exp(arb_t x, const fmpz_t y, const fmpz_t exp)
{
    arf_set_fmpz_2exp(arb_midref(x), y, exp);
    mag_zero(arb_radref(x));
}

void arb_set_round_fmpz_2exp(arb_t y, const fmpz_t x, const fmpz_t exp, long prec);

void arb_set_round_fmpz(arb_t y, const fmpz_t x, long prec);

ARB_INLINE int
arb_is_one(const arb_t f)
{
    return arf_is_one(arb_midref(f)) && mag_is_zero(arb_radref(f));
}

ARB_INLINE void
arb_one(arb_t f)
{
    arb_set_ui(f, 1UL);
}

ARB_INLINE void
arb_print(const arb_t x)
{
    arf_print(arb_midref(x));
    printf(" +/- ");
    mag_print(arb_radref(x));
}

ARB_INLINE void
arb_printd(const arb_t x, long digits)
{
    arf_printd(arb_midref(x), FLINT_MAX(digits, 1));
    printf(" +/- ");
    mag_printd(arb_radref(x), 5);
}

ARB_INLINE void
arb_mul_2exp_si(arb_t y, const arb_t x, long e)
{
    arf_mul_2exp_si(arb_midref(y), arb_midref(x), e);
    mag_mul_2exp_si(arb_radref(y), arb_radref(x), e);
}

ARB_INLINE void
arb_mul_2exp_fmpz(arb_t y, const arb_t x, const fmpz_t e)
{
    arf_mul_2exp_fmpz(arb_midref(y), arb_midref(x), e);
    mag_mul_2exp_fmpz(arb_radref(y), arb_radref(x), e);
}

ARB_INLINE int
arb_is_int(const arb_t x)
{
    return mag_is_zero(arb_radref(x)) &&
           arf_is_int(arb_midref(x));
}

ARB_INLINE int
arb_contains_zero(const arb_t x)
{
    return arf_cmpabs_mag(arb_midref(x), arb_radref(x)) <= 0;
}

ARB_INLINE int
arb_is_nonzero(const arb_t x)
{
    return !arb_contains_zero(x);
}

ARB_INLINE int
arb_is_positive(const arb_t x)
{
    return (arf_sgn(arb_midref(x)) > 0) &&
        (arf_mag_cmpabs(arb_radref(x), arb_midref(x)) < 0) &&
         !arf_is_nan(arb_midref(x));
}

ARB_INLINE int
arb_is_nonnegative(const arb_t x)
{
    return (arf_sgn(arb_midref(x)) >= 0) &&
        (arf_mag_cmpabs(arb_radref(x), arb_midref(x)) <= 0) &&
         !arf_is_nan(arb_midref(x));
}

ARB_INLINE int
arb_is_negative(const arb_t x)
{
    return (arf_sgn(arb_midref(x)) < 0) &&
        (arf_mag_cmpabs(arb_radref(x), arb_midref(x)) < 0) &&
         !arf_is_nan(arb_midref(x));
}

ARB_INLINE int
arb_is_nonpositive(const arb_t x)
{
    return (arf_sgn(arb_midref(x)) <= 0) &&
        (arf_mag_cmpabs(arb_radref(x), arb_midref(x)) <= 0) &&
         !arf_is_nan(arb_midref(x));
}

ARB_INLINE int
arb_contains_negative(const arb_t x)
{
    return (arf_sgn(arb_midref(x)) < 0) ||
        (arf_mag_cmpabs(arb_radref(x), arb_midref(x)) > 0)
        || arf_is_nan(arb_midref(x));
}

ARB_INLINE int
arb_contains_nonpositive(const arb_t x)
{
    return (arf_sgn(arb_midref(x)) <= 0) ||
        (arf_mag_cmpabs(arb_radref(x), arb_midref(x)) >= 0)
        || arf_is_nan(arb_midref(x));
}

ARB_INLINE int
arb_contains_positive(const arb_t x)
{
    return (arf_sgn(arb_midref(x)) > 0) ||
        (arf_mag_cmpabs(arb_radref(x), arb_midref(x)) > 0)
        || arf_is_nan(arb_midref(x));
}

ARB_INLINE int
arb_contains_nonnegative(const arb_t x)
{
    return (arf_sgn(arb_midref(x)) >= 0) ||
        (arf_mag_cmpabs(arb_radref(x), arb_midref(x)) >= 0)
        || arf_is_nan(arb_midref(x));
}

void arb_get_mag_lower(mag_t z, const arb_t x);

void arb_get_mag_lower_nonnegative(mag_t z, const arb_t x);

ARB_INLINE void
arb_get_mag(mag_t z, const arb_t x)
{
    mag_t t;
    mag_init_set_arf(t, arb_midref(x));
    mag_add(z, t, arb_radref(x));
    mag_clear(t);
}

ARB_INLINE void
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

ARB_INLINE void
arb_get_abs_lbound_arf(arf_t u, const arb_t x, long prec)
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

long arb_rel_error_bits(const arb_t x);

ARB_INLINE long
arb_rel_accuracy_bits(const arb_t x)
{
    return -arb_rel_error_bits(x);
}

ARB_INLINE long
arb_bits(const arb_t x)
{
    return arf_bits(arb_midref(x));
}

void arb_randtest_exact(arb_t x, flint_rand_t state, long prec, long mag_bits);

void arb_randtest_wide(arb_t x, flint_rand_t state, long prec, long mag_bits);

void arb_randtest_precise(arb_t x, flint_rand_t state, long prec, long mag_bits);

void arb_randtest(arb_t x, flint_rand_t state, long prec, long mag_bits);

void arb_randtest_special(arb_t x, flint_rand_t state, long prec, long mag_bits);

void arb_add_error_arf(arb_t x, const arf_t err);

void arb_add_error_2exp_si(arb_t x, long err);

void arb_add_error_2exp_fmpz(arb_t x, const fmpz_t err);

void arb_add_error(arb_t x, const arb_t error);

/* TODO: document */
ARB_INLINE void
arb_add_error_mag(arb_t x, const mag_t err)
{
    mag_add(arb_radref(x), arb_radref(x), err);
}

int arb_contains_arf(const arb_t x, const arf_t y);

int arb_contains_fmpq(const arb_t x, const fmpq_t y);

int arb_contains_fmpz(const arb_t x, const fmpz_t y);

int arb_contains_si(const arb_t x, long y);

int arb_contains_mpfr(const arb_t x, const mpfr_t y);

int arb_overlaps(const arb_t x, const arb_t y);

int arb_contains(const arb_t x, const arb_t y);

void arb_get_interval_fmpz_2exp(fmpz_t a, fmpz_t b, fmpz_t exp, const arb_t x);
int arb_get_unique_fmpz(fmpz_t z, const arb_t x);

void arb_floor(arb_t z, const arb_t x, long prec);
void arb_ceil(arb_t z, const arb_t x, long prec);

void arb_set_interval_arf(arb_t x, const arf_t a, const arf_t b, long prec);
void arb_set_interval_mpfr(arb_t x, const mpfr_t a, const mpfr_t b, long prec);

void arb_get_interval_arf(arf_t a, arf_t b, const arb_t x, long prec);
void arb_get_interval_mpfr(mpfr_t a, mpfr_t b, const arb_t x);

void arb_union(arb_t z, const arb_t x, const arb_t y, long prec);
void arb_get_rand_fmpq(fmpq_t q, flint_rand_t state, const arb_t x, long bits);

void arb_add(arb_t z, const arb_t x, const arb_t y, long prec);
void arb_add_arf(arb_t z, const arb_t x, const arf_t y, long prec);
void arb_add_ui(arb_t z, const arb_t x, ulong y, long prec);
void arb_add_si(arb_t z, const arb_t x, long y, long prec);
void arb_add_fmpz(arb_t z, const arb_t x, const fmpz_t y, long prec);
void arb_add_fmpz_2exp(arb_t z, const arb_t x, const fmpz_t man, const fmpz_t exp, long prec);

void arb_sub(arb_t z, const arb_t x, const arb_t y, long prec);
void arb_sub_arf(arb_t z, const arb_t x, const arf_t y, long prec);
void arb_sub_ui(arb_t z, const arb_t x, ulong y, long prec);
void arb_sub_si(arb_t z, const arb_t x, long y, long prec);
void arb_sub_fmpz(arb_t z, const arb_t x, const fmpz_t y, long prec);

void arb_mul(arb_t z, const arb_t x, const arb_t y, long prec);
void arb_mul_arf(arb_t z, const arb_t x, const arf_t y, long prec);
void arb_mul_si(arb_t z, const arb_t x, long y, long prec);
void arb_mul_ui(arb_t z, const arb_t x, ulong y, long prec);
void arb_mul_fmpz(arb_t z, const arb_t x, const fmpz_t y, long prec);

void arb_addmul(arb_t z, const arb_t x, const arb_t y, long prec);
void arb_addmul_arf(arb_t z, const arb_t x, const arf_t y, long prec);
void arb_addmul_si(arb_t z, const arb_t x, long y, long prec);
void arb_addmul_ui(arb_t z, const arb_t x, ulong y, long prec);
void arb_addmul_fmpz(arb_t z, const arb_t x, const fmpz_t y, long prec);

void arb_submul(arb_t z, const arb_t x, const arb_t y, long prec);
void arb_submul_arf(arb_t z, const arb_t x, const arf_t y, long prec);
void arb_submul_si(arb_t z, const arb_t x, long y, long prec);
void arb_submul_ui(arb_t z, const arb_t x, ulong y, long prec);
void arb_submul_fmpz(arb_t z, const arb_t x, const fmpz_t y, long prec);

void arb_div(arb_t z, const arb_t x, const arb_t y, long prec);
void arb_div_arf(arb_t z, const arb_t x, const arf_t y, long prec);
void arb_div_si(arb_t z, const arb_t x, long y, long prec);
void arb_div_ui(arb_t z, const arb_t x, ulong y, long prec);
void arb_div_fmpz(arb_t z, const arb_t x, const fmpz_t y, long prec);
void arb_fmpz_div_fmpz(arb_t z, const fmpz_t x, const fmpz_t y, long prec);
void arb_ui_div(arb_t z, ulong x, const arb_t y, long prec);

ARB_INLINE void
arb_inv(arb_t y, const arb_t x, long prec)
{
    arb_ui_div(y, 1, x, prec);
}

ARB_INLINE void
arb_set_fmpq(arb_t y, const fmpq_t x, long prec)
{
    arb_fmpz_div_fmpz(y, fmpq_numref(x), fmpq_denref(x), prec);
}

void arb_sqrt(arb_t z, const arb_t x, long prec);
void arb_sqrt_arf(arb_t z, const arf_t x, long prec);
void arb_sqrt_fmpz(arb_t z, const fmpz_t x, long prec);
void arb_sqrt_ui(arb_t z, ulong x, long prec);

void arb_sqrtpos(arb_t z, const arb_t x, long prec);
void arb_hypot(arb_t z, const arb_t x, const arb_t y, long prec);

void arb_rsqrt(arb_t z, const arb_t x, long prec);
void arb_rsqrt_ui(arb_t z, ulong x, long prec);

void arb_pow_fmpz_binexp(arb_t y, const arb_t b, const fmpz_t e, long prec);
void arb_pow_fmpz(arb_t y, const arb_t b, const fmpz_t e, long prec);
void arb_pow_ui(arb_t y, const arb_t b, ulong e, long prec);
void arb_ui_pow_ui(arb_t y, ulong b, ulong e, long prec);
void arb_si_pow_ui(arb_t y, long b, ulong e, long prec);
void arb_pow_fmpq(arb_t y, const arb_t x, const fmpq_t a, long prec);

void arb_div_2expm1_ui(arb_t z, const arb_t x, ulong n, long prec);
void arb_pow(arb_t z, const arb_t x, const arb_t y, long prec);
void arb_root(arb_t z, const arb_t x, ulong k, long prec);
void arb_log(arb_t z, const arb_t x, long prec);
void arb_log_arf(arb_t z, const arf_t x, long prec);
void arb_log_ui(arb_t z, ulong x, long prec);
void arb_log_fmpz(arb_t z, const fmpz_t x, long prec);
void arb_exp(arb_t z, const arb_t x, long prec);
void arb_expm1(arb_t z, const arb_t x, long prec);
void arb_sin(arb_t s, const arb_t x, long prec);
void arb_cos(arb_t c, const arb_t x, long prec);
void arb_sin_cos(arb_t s, arb_t c, const arb_t x, long prec);
void arb_sin_pi(arb_t s, const arb_t x, long prec);
void arb_cos_pi(arb_t c, const arb_t x, long prec);
void arb_sin_cos_pi(arb_t s, arb_t c, const arb_t x, long prec);
void arb_tan(arb_t y, const arb_t x, long prec);
void arb_cot(arb_t y, const arb_t x, long prec);
void arb_tan_pi(arb_t y, const arb_t x, long prec);
void arb_cot_pi(arb_t y, const arb_t x, long prec);
void _arb_sin_pi_fmpq_algebraic(arb_t s, ulong p, ulong q, long prec);
void _arb_cos_pi_fmpq_algebraic(arb_t c, ulong p, ulong q, long prec);
void _arb_sin_cos_pi_fmpq_algebraic(arb_t s, arb_t c, ulong p, ulong q, long prec);
void arb_sin_cos_pi_fmpq(arb_t s, arb_t c, const fmpq_t x, long prec);
void arb_sin_pi_fmpq(arb_t s, const fmpq_t x, long prec);
void arb_cos_pi_fmpq(arb_t c, const fmpq_t x, long prec);
void arb_sinh(arb_t z, const arb_t x, long prec);
void arb_cosh(arb_t z, const arb_t x, long prec);
void arb_sinh_cosh(arb_t s, arb_t c, const arb_t x, long prec);
void arb_tanh(arb_t y, const arb_t x, long prec);
void arb_coth(arb_t y, const arb_t x, long prec);
void arb_atan_arf(arb_t z, const arf_t x, long prec);
void arb_atan(arb_t z, const arb_t x, long prec);
void arb_atan2(arb_t z, const arb_t b, const arb_t a, long prec);
void arb_asin(arb_t z, const arb_t x, long prec);
void arb_acos(arb_t z, const arb_t x, long prec);
void arb_fac_ui(arb_t z, ulong n, long prec);
void arb_bin_ui(arb_t z, const arb_t n, ulong k, long prec);
void arb_bin_uiui(arb_t z, ulong n, ulong k, long prec);
void arb_fib_fmpz(arb_t z, const fmpz_t n, long prec);
void arb_fib_ui(arb_t z, ulong n, long prec);
void arb_const_pi(arb_t z, long prec);
void arb_const_sqrt_pi(arb_t z, long prec);
void arb_const_log_sqrt2pi(arb_t z, long prec);
void arb_const_log2(arb_t z, long prec);
void arb_const_log10(arb_t z, long prec);
void arb_const_euler(arb_t z, long prec);
void arb_const_catalan(arb_t z, long prec);
void arb_const_e(arb_t z, long prec);
void arb_const_khinchin(arb_t z, long prec);
void arb_const_glaisher(arb_t z, long prec);
void arb_agm(arb_t z, const arb_t x, const arb_t y, long prec);
void arb_lgamma(arb_t z, const arb_t x, long prec);
void arb_rgamma(arb_t z, const arb_t x, long prec);
void arb_gamma(arb_t z, const arb_t x, long prec);
void arb_gamma_fmpq(arb_t z, const fmpq_t x, long prec);
void arb_gamma_fmpz(arb_t z, const fmpz_t x, long prec);
void arb_digamma(arb_t y, const arb_t x, long prec);
void arb_zeta(arb_t z, const arb_t s, long prec);
void arb_zeta_ui(arb_t z, ulong n, long prec);
void arb_bernoulli_ui(arb_t z, ulong n, long prec);

void arb_rising_ui_bs(arb_t y, const arb_t x, ulong n, long prec);
void arb_rising_ui_rs(arb_t y, const arb_t x, ulong n, ulong m, long prec);
void arb_rising_ui_rec(arb_t y, const arb_t x, ulong n, long prec);
void arb_rising_ui(arb_t z, const arb_t x, ulong n, long prec);
void arb_rising_fmpq_ui(arb_t y, const fmpq_t x, ulong n, long prec);

void arb_rising2_ui_rs(arb_t u, arb_t v, const arb_t x, ulong n, ulong m, long prec);
void arb_rising2_ui_bs(arb_t u, arb_t v, const arb_t x, ulong n, long prec);
void arb_rising2_ui(arb_t u, arb_t v, const arb_t x, ulong n, long prec);

void arb_log_ui_from_prev(arb_t s, ulong k, arb_t log_prev, ulong prev, long prec);

void arb_const_apery(arb_t s, long prec);

void arb_zeta_ui_asymp(arb_t x, ulong s, long prec);
void arb_zeta_ui_borwein_bsplit(arb_t x, ulong s, long prec);
void arb_zeta_ui_euler_product(arb_t z, ulong s, long prec);
void arb_zeta_ui_bernoulli(arb_t x, ulong n, long prec);
void arb_zeta_ui_vec_borwein(arb_ptr z, ulong start, long num, ulong step, long prec);
void arb_zeta_ui(arb_t x, ulong n, long prec);
void arb_zeta_ui_vec_even(arb_ptr x, ulong start, long num, long prec);
void arb_zeta_ui_vec_odd(arb_ptr x, ulong start, long num, long prec);
void arb_zeta_ui_vec(arb_ptr x, ulong start, long num, long prec);
void arb_bernoulli_ui(arb_t b, ulong n, long prec);
void arb_bernoulli_ui_zeta(arb_t b, ulong n, long prec);

void arb_polylog(arb_t w, const arb_t s, const arb_t z, long prec);
void arb_polylog_si(arb_t w, long s, const arb_t z, long prec);

void arb_chebyshev_t_ui(arb_t a, ulong n, const arb_t x, long prec);
void arb_chebyshev_t2_ui(arb_t a, arb_t b, ulong n, const arb_t x, long prec);
void arb_chebyshev_u_ui(arb_t a, ulong n, const arb_t x, long prec);
void arb_chebyshev_u2_ui(arb_t a, arb_t b, ulong n, const arb_t x, long prec);

#define ARB_DEF_CACHED_CONSTANT(name, comp_func) \
    TLS_PREFIX long name ## _cached_prec = 0; \
    TLS_PREFIX arb_t name ## _cached_value; \
    void name ## _cleanup(void) \
    { \
        arb_clear(name ## _cached_value); \
        name ## _cached_prec = 0; \
    } \
    void name(arb_t x, long prec) \
    { \
        if (name ## _cached_prec < prec) \
        { \
            if (name ## _cached_prec == 0) \
            { \
                arb_init(name ## _cached_value); \
                flint_register_cleanup_function(name ## _cleanup); \
            } \
            comp_func(name ## _cached_value, prec); \
            name ## _cached_prec = prec; \
        } \
        arb_set_round(x, name ## _cached_value, prec); \
    }

/* vector functions */

ARB_INLINE void
_arb_vec_zero(arb_ptr A, long n)
{
    long i;
    for (i = 0; i < n; i++)
        arb_zero(A + i);
}

ARB_INLINE int
_arb_vec_is_zero(arb_srcptr vec, long len)
{
    long i;
    for (i = 0; i < len; i++)
        if (!arb_is_zero(vec + i))
            return 0;
    return 1;
}

ARB_INLINE int
_arb_vec_is_finite(arb_srcptr x, long len)
{
    long i;

    for (i = 0; i < len; i++)
        if (!arb_is_finite(x + i))
            return 0;

    return 1;
}

ARB_INLINE void
_arb_vec_set(arb_ptr res, arb_srcptr vec, long len)
{
    long i;
    for (i = 0; i < len; i++)
        arb_set(res + i, vec + i);
}

ARB_INLINE void
_arb_vec_set_round(arb_ptr res, arb_srcptr vec, long len, long prec)
{
    long i;
    for (i = 0; i < len; i++)
        arb_set_round(res + i, vec + i, prec);
}

ARB_INLINE void
_arb_vec_swap(arb_ptr res, arb_ptr vec, long len)
{
    long i;
    for (i = 0; i < len; i++)
        arb_swap(res + i, vec + i);
}

ARB_INLINE void
_arb_vec_neg(arb_ptr B, arb_srcptr A, long n)
{
    long i;
    for (i = 0; i < n; i++)
        arb_neg(B + i, A + i);
}

ARB_INLINE void
_arb_vec_sub(arb_ptr C, arb_srcptr A,
    arb_srcptr B, long n, long prec)
{
    long i;
    for (i = 0; i < n; i++)
        arb_sub(C + i, A + i, B + i, prec);
}

ARB_INLINE void
_arb_vec_add(arb_ptr C, arb_srcptr A,
    arb_srcptr B, long n, long prec)
{
    long i;
    for (i = 0; i < n; i++)
        arb_add(C + i, A + i, B + i, prec);
}

ARB_INLINE void
_arb_vec_scalar_mul(arb_ptr res, arb_srcptr vec,
    long len, const arb_t c, long prec)
{
    long i;
    for (i = 0; i < len; i++)
        arb_mul(res + i, vec + i, c, prec);
}

ARB_INLINE void
_arb_vec_scalar_div(arb_ptr res, arb_srcptr vec,
    long len, const arb_t c, long prec)
{
    long i;
    for (i = 0; i < len; i++)
        arb_div(res + i, vec + i, c, prec);
}

ARB_INLINE void
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

ARB_INLINE void
_arb_vec_scalar_mul_2exp_si(arb_ptr res, arb_srcptr src, long len, long c)
{
    long i;
    for (i = 0; i < len; i++)
        arb_mul_2exp_si(res + i, src + i, c);
}

ARB_INLINE void
_arb_vec_scalar_addmul(arb_ptr res, arb_srcptr vec,
    long len, const arb_t c, long prec)
{
    long i;
    for (i = 0; i < len; i++)
        arb_addmul(res + i, vec + i, c, prec);
}

ARB_INLINE void
_arb_vec_dot(arb_t res, arb_srcptr vec1, arb_srcptr vec2, long len2, long prec)
{
    long i;
    arb_zero(res);
    for (i = 0; i < len2; i++)
        arb_addmul(res, vec1 + i, vec2 + i, prec);
}

ARB_INLINE void
_arb_vec_norm(arb_t res, arb_srcptr vec, long len, long prec)
{
    long i;
    arb_zero(res);
    for (i = 0; i < len; i++)
        arb_addmul(res, vec + i, vec + i, prec);
}

/* TODO: mag version? */
ARB_INLINE void
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

ARB_INLINE long
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

ARB_INLINE void
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

ARB_INLINE void
_arb_vec_add_error_arf_vec(arb_ptr res, arf_srcptr err, long len)
{
    long i;
    for (i = 0; i < len; i++)
        arb_add_error_arf(res + i, err + i);
}

ARB_INLINE void
_arb_vec_add_error_mag_vec(arb_ptr res, mag_srcptr err, long len)
{
    long i;
    for (i = 0; i < len; i++)
        mag_add(arb_radref(res + i), arb_radref(res + i), err + i);
}

ARB_INLINE void
_arb_vec_indeterminate(arb_ptr vec, long len)
{
    long i;
    for (i = 0; i < len; i++)
        arb_indeterminate(vec + i);
}

ARB_INLINE void
_arb_vec_trim(arb_ptr res, arb_srcptr vec, long len)
{
    long i;
    for (i = 0; i < len; i++)
        arb_trim(res + i, vec + i);
}

/* arctangent implementation */

#define ARB_ATAN_TAB1_BITS 8
#define ARB_ATAN_TAB1_PREC 512
#define ARB_ATAN_TAB1_LIMBS (ARB_ATAN_TAB1_PREC / FLINT_BITS)

#define ARB_ATAN_TAB21_BITS 5
#define ARB_ATAN_TAB22_BITS 5
#define ARB_ATAN_TAB2_PREC 4608
#define ARB_ATAN_TAB2_LIMBS (ARB_ATAN_TAB2_PREC / FLINT_BITS)

extern const mp_limb_t arb_atan_tab1[1 << ARB_ATAN_TAB1_BITS][ARB_ATAN_TAB1_LIMBS];
extern const mp_limb_t arb_atan_tab21[1 << ARB_ATAN_TAB21_BITS][ARB_ATAN_TAB2_LIMBS];
extern const mp_limb_t arb_atan_tab22[1 << ARB_ATAN_TAB22_BITS][ARB_ATAN_TAB2_LIMBS];
extern const mp_limb_t arb_atan_pi2_minus_one[ARB_ATAN_TAB2_LIMBS];

void
_arb_atan_taylor_naive(mp_ptr y, mp_limb_t * error,
    mp_srcptr x, mp_size_t xn, ulong N, int alternating);

void _arb_atan_taylor_rs(mp_ptr y, mp_limb_t * error,
    mp_srcptr x, mp_size_t xn, ulong N, int alternating);

/* logarithm implementation */

#define ARB_LOG_TAB11_BITS 7
#define ARB_LOG_TAB12_BITS 7
#define ARB_LOG_TAB1_PREC 512
#define ARB_LOG_TAB1_LIMBS (ARB_LOG_TAB1_PREC / FLINT_BITS)

#define ARB_LOG_TAB21_BITS 5
#define ARB_LOG_TAB22_BITS 5
#define ARB_LOG_TAB2_PREC 4608
#define ARB_LOG_TAB2_LIMBS (ARB_LOG_TAB2_PREC / FLINT_BITS)

extern const mp_limb_t arb_log_tab11[1 << ARB_LOG_TAB11_BITS][ARB_LOG_TAB1_LIMBS];
extern const mp_limb_t arb_log_tab12[1 << ARB_LOG_TAB12_BITS][ARB_LOG_TAB1_LIMBS];
extern const mp_limb_t arb_log_tab21[1 << ARB_LOG_TAB21_BITS][ARB_LOG_TAB2_LIMBS];
extern const mp_limb_t arb_log_tab22[1 << ARB_LOG_TAB22_BITS][ARB_LOG_TAB2_LIMBS];
extern const mp_limb_t arb_log_log2_tab[ARB_LOG_TAB2_LIMBS];

/* exponential implementation */

/* only goes up to log(2) * 256 */
#define ARB_EXP_TAB1_NUM 178
#define ARB_EXP_TAB1_BITS 8
#define ARB_EXP_TAB1_PREC 512
#define ARB_EXP_TAB1_LIMBS (ARB_EXP_TAB1_PREC / FLINT_BITS)

/* only goes up to log(2) * 32 */
#define ARB_EXP_TAB21_NUM 23
#define ARB_EXP_TAB21_BITS 5
#define ARB_EXP_TAB22_NUM (1 << ARB_EXP_TAB22_BITS)
#define ARB_EXP_TAB22_BITS 5
#define ARB_EXP_TAB2_PREC 4608
#define ARB_EXP_TAB2_LIMBS (ARB_EXP_TAB2_PREC / FLINT_BITS)

extern const mp_limb_t arb_exp_tab1[ARB_EXP_TAB1_NUM][ARB_EXP_TAB1_LIMBS];
extern const mp_limb_t arb_exp_tab21[ARB_EXP_TAB21_NUM][ARB_EXP_TAB2_LIMBS];
extern const mp_limb_t arb_exp_tab22[ARB_EXP_TAB22_NUM][ARB_EXP_TAB2_LIMBS];

void _arb_exp_taylor_naive(mp_ptr y, mp_limb_t * error,
    mp_srcptr x, mp_size_t xn, ulong N);

void _arb_exp_taylor_rs(mp_ptr y, mp_limb_t * error,
    mp_srcptr x, mp_size_t xn, ulong N);

void arb_exp_arf_bb(arb_t z, const arf_t x, long prec, int minus_one);

int _arb_get_mpn_fixed_mod_log2(mp_ptr w, fmpz_t q, mp_limb_t * error,
                                                const arf_t x, mp_size_t wn);

long _arb_exp_taylor_bound(long mag, long prec);

void _arb_exp_sum_bs_powtab(fmpz_t T, fmpz_t Q, mp_bitcnt_t * Qexp,
    const fmpz_t x, mp_bitcnt_t r, long N);

void _arb_exp_sum_bs_simple(fmpz_t T, fmpz_t Q, mp_bitcnt_t * Qexp,
    const fmpz_t x, mp_bitcnt_t r, long N);

/* sin/cos implementation */

/* only goes up to (pi/4) * 256 */
#define ARB_SIN_COS_TAB1_NUM 203
#define ARB_SIN_COS_TAB1_BITS 8
#define ARB_SIN_COS_TAB1_PREC 512
#define ARB_SIN_COS_TAB1_LIMBS (ARB_SIN_COS_TAB1_PREC / FLINT_BITS)

/* only goes up to (pi/4) * 32 */
#define ARB_SIN_COS_TAB21_NUM 26
#define ARB_SIN_COS_TAB21_BITS 5
#define ARB_SIN_COS_TAB22_NUM (1 << ARB_SIN_COS_TAB22_BITS)
#define ARB_SIN_COS_TAB22_BITS 5
#define ARB_SIN_COS_TAB2_PREC 4608
#define ARB_SIN_COS_TAB2_LIMBS (ARB_SIN_COS_TAB2_PREC / FLINT_BITS)

extern const mp_limb_t arb_sin_cos_tab1[2 * ARB_SIN_COS_TAB1_NUM][ARB_SIN_COS_TAB1_LIMBS];
extern const mp_limb_t arb_sin_cos_tab21[2 * ARB_SIN_COS_TAB21_NUM][ARB_SIN_COS_TAB2_LIMBS];
extern const mp_limb_t arb_sin_cos_tab22[2 * ARB_SIN_COS_TAB22_NUM][ARB_SIN_COS_TAB2_LIMBS];

#define ARB_PI4_TAB_LIMBS (4608 / FLINT_BITS)
extern const mp_limb_t arb_pi4_tab[ARB_PI4_TAB_LIMBS];

void _arb_sin_cos_taylor_naive(mp_ptr ysin, mp_ptr ycos, mp_limb_t * error,
    mp_srcptr x, mp_size_t xn, ulong N);

void _arb_sin_cos_taylor_rs(mp_ptr ysin, mp_ptr ycos,
    mp_limb_t * error, mp_srcptr x, mp_size_t xn, ulong N,
    int sinonly, int alternating);

int _arb_get_mpn_fixed_mod_pi4(mp_ptr w, fmpz_t q, int * octant,
    mp_limb_t * error, const arf_t x, mp_size_t wn);

ARB_INLINE mp_bitcnt_t
_arb_mpn_leading_zeros(mp_srcptr d, mp_size_t n)
{
    mp_limb_t t;
    mp_size_t zero_limbs;
    mp_bitcnt_t bits;

    zero_limbs = 0;

    while (1)
    {
        t = d[n - zero_limbs - 1];

        if (t != 0)
        {
            count_leading_zeros(bits, t);
            return bits + FLINT_BITS * zero_limbs;
        }

        zero_limbs++;

        if (zero_limbs == n)
            return FLINT_BITS * zero_limbs;
    }
}

void _arb_atan_sum_bs_simple(fmpz_t T, fmpz_t Q, mp_bitcnt_t * Qexp,
    const fmpz_t x, mp_bitcnt_t r, long N);

void _arb_atan_sum_bs_powtab(fmpz_t T, fmpz_t Q, mp_bitcnt_t * Qexp,
    const fmpz_t x, mp_bitcnt_t r, long N);

void arb_atan_arf_bb(arb_t z, const arf_t x, long prec);

#ifdef __cplusplus
}
#endif

#endif

