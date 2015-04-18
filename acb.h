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

#ifndef ACB_H
#define ACB_H

#ifdef ACB_INLINES_C
#define ACB_INLINE
#else
#define ACB_INLINE static __inline__
#endif

#include "arf.h"
#include "arb.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    arb_struct real;
    arb_struct imag;
}
acb_struct;

typedef acb_struct acb_t[1];
typedef acb_struct * acb_ptr;
typedef const acb_struct * acb_srcptr;

#define acb_realref(x) (&(x)->real)
#define acb_imagref(x) (&(x)->imag)


ACB_INLINE void
acb_init(acb_t x)
{
    arb_init(acb_realref(x));
    arb_init(acb_imagref(x));
}

ACB_INLINE void
acb_clear(acb_t x)
{
    arb_clear(acb_realref(x));
    arb_clear(acb_imagref(x));
}

ACB_INLINE acb_ptr
_acb_vec_init(long n)
{
    long i;
    acb_ptr v = (acb_ptr) flint_malloc(sizeof(acb_struct) * n);

    for (i = 0; i < n; i++)
        acb_init(v + i);

    return v;
}

ACB_INLINE void
_acb_vec_clear(acb_ptr v, long n)
{
    long i;
    for (i = 0; i < n; i++)
        acb_clear(v + i);
    flint_free(v);
}

ACB_INLINE int
acb_is_zero(const acb_t z)
{
    return arb_is_zero(acb_realref(z)) && arb_is_zero(acb_imagref(z));
}

ACB_INLINE int
acb_is_one(const acb_t z)
{
    return arb_is_one(acb_realref(z)) && arb_is_zero(acb_imagref(z));
}

ACB_INLINE int
acb_is_exact(const acb_t z)
{
    return arb_is_exact(acb_realref(z)) && arb_is_exact(acb_imagref(z));
}

ACB_INLINE int
acb_is_int(const acb_t z)
{
    return arb_is_zero(acb_imagref(z)) && arb_is_int(acb_realref(z));
}

ACB_INLINE void
acb_zero(acb_t z)
{
    arb_zero(acb_realref(z));
    arb_zero(acb_imagref(z));
}

ACB_INLINE void
acb_one(acb_t z)
{
    arb_one(acb_realref(z));
    arb_zero(acb_imagref(z));
}

ACB_INLINE void
acb_onei(acb_t z)
{
    arb_zero(acb_realref(z));
    arb_one(acb_imagref(z));
}

ACB_INLINE void
acb_set(acb_t z, const acb_t x)
{
    arb_set(acb_realref(z), acb_realref(x));
    arb_set(acb_imagref(z), acb_imagref(x));
}

ACB_INLINE void
acb_set_round(acb_t z, const acb_t x, long prec)
{
    arb_set_round(acb_realref(z), acb_realref(x), prec);
    arb_set_round(acb_imagref(z), acb_imagref(x), prec);
}

ACB_INLINE void
acb_neg_round(acb_t z, const acb_t x, long prec)
{
    arb_neg_round(acb_realref(z), acb_realref(x), prec);
    arb_neg_round(acb_imagref(z), acb_imagref(x), prec);
}

ACB_INLINE void
acb_swap(acb_t z, acb_t x)
{
    arb_swap(acb_realref(z), acb_realref(x));
    arb_swap(acb_imagref(z), acb_imagref(x));
}

ACB_INLINE int
acb_equal(const acb_t x, const acb_t y)
{
    return arb_equal(acb_realref(x), acb_realref(y)) &&
           arb_equal(acb_imagref(x), acb_imagref(y));
}

ACB_INLINE int
acb_eq(const acb_t x, const acb_t y)
{
    return arb_eq(acb_realref(x), acb_realref(y)) &&
           arb_eq(acb_imagref(x), acb_imagref(y));
}

ACB_INLINE int
acb_ne(const acb_t x, const acb_t y)
{
    return arb_ne(acb_realref(x), acb_realref(y)) ||
           arb_ne(acb_imagref(x), acb_imagref(y));
}

ACB_INLINE int
acb_overlaps(const acb_t x, const acb_t y)
{
    return arb_overlaps(acb_realref(x), acb_realref(y)) &&
            arb_overlaps(acb_imagref(x), acb_imagref(y));
}

ACB_INLINE int
acb_contains_zero(const acb_t x)
{
    return arb_contains_zero(acb_realref(x)) &&
            arb_contains_zero(acb_imagref(x));
}

ACB_INLINE int
acb_contains_fmpq(const acb_t x, const fmpq_t y)
{
    return arb_contains_fmpq(acb_realref(x), y) &&
            arb_contains_zero(acb_imagref(x));
}

ACB_INLINE int
acb_contains_fmpz(const acb_t x, const fmpz_t y)
{
    return arb_contains_fmpz(acb_realref(x), y) &&
            arb_contains_zero(acb_imagref(x));
}

ACB_INLINE int
acb_contains(const acb_t x, const acb_t y)
{
    return arb_contains(acb_realref(x), acb_realref(y)) &&
            arb_contains(acb_imagref(x), acb_imagref(y));
}

ACB_INLINE void
acb_set_ui(acb_t z, ulong c)
{
    arb_set_ui(acb_realref(z), c);
    arb_zero(acb_imagref(z));
}

ACB_INLINE void
acb_set_si(acb_t z, long c)
{
    arb_set_si(acb_realref(z), c);
    arb_zero(acb_imagref(z));
}

ACB_INLINE void
acb_set_fmpz(acb_t z, const fmpz_t c)
{
    arb_set_fmpz(acb_realref(z), c);
    arb_zero(acb_imagref(z));
}

ACB_INLINE void
acb_set_round_fmpz(acb_t z, const fmpz_t y, long prec)
{
    arb_set_round_fmpz(acb_realref(z), y, prec);
    arb_zero(acb_imagref(z));
}

int acb_get_unique_fmpz(fmpz_t z, const acb_t x);

ACB_INLINE void
acb_set_fmpq(acb_t z, const fmpq_t c, long prec)
{
    arb_set_fmpq(acb_realref(z), c, prec);
    arb_zero(acb_imagref(z));
}

ACB_INLINE void
acb_set_arb(acb_t z, const arb_t c)
{
    arb_set(acb_realref(z), c);
    arb_zero(acb_imagref(z));
}

ACB_INLINE void
acb_set_round_arb(acb_t z, const arb_t x, long prec)
{
    arb_set_round(acb_realref(z), x, prec);
    arb_zero(acb_imagref(z));
}

ACB_INLINE void
acb_trim(acb_t z, const acb_t x)
{
    arb_trim(acb_realref(z), acb_realref(x));
    arb_trim(acb_imagref(z), acb_imagref(x));
}

ACB_INLINE void
acb_add_error_arf(acb_t x, const arf_t err)
{
    arb_add_error_arf(acb_realref(x), err);
    arb_add_error_arf(acb_imagref(x), err);
}

/* TODO: document */
ACB_INLINE void
acb_add_error_mag(acb_t x, const mag_t err)
{
    arb_add_error_mag(acb_realref(x), err);
    arb_add_error_mag(acb_imagref(x), err);
}

void acb_get_mag(mag_t z, const acb_t x);

void acb_get_mag_lower(mag_t z, const acb_t x);

ACB_INLINE void
acb_get_abs_ubound_arf(arf_t u, const acb_t z, long prec)
{
    if (arb_is_zero(acb_imagref(z)))
    {
        arb_get_abs_ubound_arf(u, acb_realref(z), prec);
    }
    else if (arb_is_zero(acb_realref(z)))
    {
        arb_get_abs_ubound_arf(u, acb_imagref(z), prec);
    }
    else
    {
        arf_t v;
        arf_init(v);

        arb_get_abs_ubound_arf(u, acb_realref(z), prec);
        arb_get_abs_ubound_arf(v, acb_imagref(z), prec);

        arf_mul(u, u, u, prec, ARF_RND_UP);
        arf_mul(v, v, v, prec, ARF_RND_UP);
        arf_add(u, u, v, prec, ARF_RND_UP);
        arf_sqrt(u, u, prec, ARF_RND_UP);

        arf_clear(v);
    }
}

ACB_INLINE void
acb_get_abs_lbound_arf(arf_t u, const acb_t z, long prec)
{
    if (arb_is_zero(acb_imagref(z)))
    {
        arb_get_abs_lbound_arf(u, acb_realref(z), prec);
    }
    else if (arb_is_zero(acb_realref(z)))
    {
        arb_get_abs_lbound_arf(u, acb_imagref(z), prec);
    }
    else
    {
        arf_t v;
        arf_init(v);

        arb_get_abs_lbound_arf(u, acb_realref(z), prec);
        arb_get_abs_lbound_arf(v, acb_imagref(z), prec);

        arf_mul(u, u, u, prec, ARF_RND_DOWN);
        arf_mul(v, v, v, prec, ARF_RND_DOWN);
        arf_add(u, u, v, prec, ARF_RND_DOWN);
        arf_sqrt(u, u, prec, ARF_RND_DOWN);

        arf_clear(v);
    }
}

ACB_INLINE void
acb_get_rad_ubound_arf(arf_t u, const acb_t z, long prec)
{
    /* fixme: this bound is very sloppy */

    if (mag_cmp(arb_radref(acb_realref(z)), arb_radref(acb_imagref(z))) >= 0)
        arf_set_mag(u, arb_radref(acb_realref(z)));
    else
        arf_set_mag(u, arb_radref(acb_imagref(z)));

    arf_mul_2exp_si(u, u, 1);
}

void acb_arg(arb_t r, const acb_t z, long prec);

ACB_INLINE void
acb_add(acb_t z, const acb_t x, const acb_t y, long prec)
{
    arb_add(acb_realref(z), acb_realref(x), acb_realref(y), prec);
    arb_add(acb_imagref(z), acb_imagref(x), acb_imagref(y), prec);
}

ACB_INLINE void
acb_sub(acb_t z, const acb_t x, const acb_t y, long prec)
{
    arb_sub(acb_realref(z), acb_realref(x), acb_realref(y), prec);
    arb_sub(acb_imagref(z), acb_imagref(x), acb_imagref(y), prec);
}

ACB_INLINE void
acb_add_ui(acb_t z, const acb_t x, ulong c, long prec)
{
    arb_add_ui(acb_realref(z), acb_realref(x), c, prec);
    arb_set_round(acb_imagref(z), acb_imagref(x), prec);
}

ACB_INLINE void
acb_sub_ui(acb_t z, const acb_t x, ulong c, long prec)
{
    arb_sub_ui(acb_realref(z), acb_realref(x), c, prec);
    arb_set_round(acb_imagref(z), acb_imagref(x), prec);
}

ACB_INLINE void
acb_add_fmpz(acb_t z, const acb_t x, const fmpz_t y, long prec)
{
    arb_add_fmpz(acb_realref(z), acb_realref(x), y, prec);
    arb_set_round(acb_imagref(z), acb_imagref(x), prec);
}

ACB_INLINE void
acb_add_arb(acb_t z, const acb_t x, const arb_t y, long prec)
{
    arb_add(acb_realref(z), acb_realref(x), y, prec);
    arb_set_round(acb_imagref(z), acb_imagref(x), prec);
}

ACB_INLINE void
acb_sub_fmpz(acb_t z, const acb_t x, const fmpz_t y, long prec)
{
    arb_sub_fmpz(acb_realref(z), acb_realref(x), y, prec);
    arb_set_round(acb_imagref(z), acb_imagref(x), prec);
}

ACB_INLINE void
acb_sub_arb(acb_t z, const acb_t x, const arb_t y, long prec)
{
    arb_sub(acb_realref(z), acb_realref(x), y, prec);
    arb_set_round(acb_imagref(z), acb_imagref(x), prec);
}

ACB_INLINE void
acb_neg(acb_t z, const acb_t x)
{
    arb_neg(acb_realref(z), acb_realref(x));
    arb_neg(acb_imagref(z), acb_imagref(x));
}

ACB_INLINE void
acb_conj(acb_t z, const acb_t x)
{
    arb_set(acb_realref(z), acb_realref(x));
    arb_neg(acb_imagref(z), acb_imagref(x));
}

ACB_INLINE void
acb_abs(arb_t u, const acb_t z, long prec)
{
    arb_hypot(u, acb_realref(z), acb_imagref(z), prec);
}

ACB_INLINE void
acb_mul_ui(acb_t z, const acb_t x, ulong y, long prec)
{
    arb_mul_ui(acb_realref(z), acb_realref(x), y, prec);
    arb_mul_ui(acb_imagref(z), acb_imagref(x), y, prec);
}

ACB_INLINE void
acb_mul_si(acb_t z, const acb_t x, long y, long prec)
{
    arb_mul_si(acb_realref(z), acb_realref(x), y, prec);
    arb_mul_si(acb_imagref(z), acb_imagref(x), y, prec);
}

ACB_INLINE void
acb_mul_fmpz(acb_t z, const acb_t x, const fmpz_t y, long prec)
{
    arb_mul_fmpz(acb_realref(z), acb_realref(x), y, prec);
    arb_mul_fmpz(acb_imagref(z), acb_imagref(x), y, prec);
}

ACB_INLINE void
acb_mul_arb(acb_t z, const acb_t x, const arb_t y, long prec)
{
    arb_mul(acb_realref(z), acb_realref(x), y, prec);
    arb_mul(acb_imagref(z), acb_imagref(x), y, prec);
}

ACB_INLINE void
acb_mul_onei(acb_t z, const acb_t x)
{
    if (z == x)
    {
        arb_swap(acb_realref(z), acb_imagref(z));
        arb_neg(acb_realref(z), acb_realref(z));
    }
    else
    {
        arb_neg(acb_realref(z), acb_imagref(x));
        arb_set(acb_imagref(z), acb_realref(x));
    }
}

void acb_mul(acb_t z, const acb_t x, const acb_t y, long prec);

void acb_mul_naive(acb_t z, const acb_t x, const acb_t y, long prec);

ACB_INLINE void
acb_mul_2exp_si(acb_t z, const acb_t x, long e)
{
    arb_mul_2exp_si(acb_realref(z), acb_realref(x), e);
    arb_mul_2exp_si(acb_imagref(z), acb_imagref(x), e);
}

ACB_INLINE void
acb_mul_2exp_fmpz(acb_t z, const acb_t x, const fmpz_t c)
{
    arb_mul_2exp_fmpz(acb_realref(z), acb_realref(x), c);
    arb_mul_2exp_fmpz(acb_imagref(z), acb_imagref(x), c);
}

void acb_addmul(acb_t z, const acb_t x, const acb_t y, long prec);

void acb_submul(acb_t z, const acb_t x, const acb_t y, long prec);

ACB_INLINE void
acb_addmul_ui(acb_t z, const acb_t x, ulong y, long prec)
{
    arb_addmul_ui(acb_realref(z), acb_realref(x), y, prec);
    arb_addmul_ui(acb_imagref(z), acb_imagref(x), y, prec);
}

ACB_INLINE void
acb_addmul_si(acb_t z, const acb_t x, long y, long prec)
{
    arb_addmul_si(acb_realref(z), acb_realref(x), y, prec);
    arb_addmul_si(acb_imagref(z), acb_imagref(x), y, prec);
}

ACB_INLINE void
acb_submul_ui(acb_t z, const acb_t x, ulong y, long prec)
{
    arb_submul_ui(acb_realref(z), acb_realref(x), y, prec);
    arb_submul_ui(acb_imagref(z), acb_imagref(x), y, prec);
}

ACB_INLINE void
acb_submul_si(acb_t z, const acb_t x, long y, long prec)
{
    arb_submul_si(acb_realref(z), acb_realref(x), y, prec);
    arb_submul_si(acb_imagref(z), acb_imagref(x), y, prec);
}

ACB_INLINE void
acb_addmul_fmpz(acb_t z, const acb_t x, const fmpz_t y, long prec)
{
    arb_addmul_fmpz(acb_realref(z), acb_realref(x), y, prec);
    arb_addmul_fmpz(acb_imagref(z), acb_imagref(x), y, prec);
}

ACB_INLINE void
acb_submul_fmpz(acb_t z, const acb_t x, const fmpz_t y, long prec)
{
    arb_submul_fmpz(acb_realref(z), acb_realref(x), y, prec);
    arb_submul_fmpz(acb_imagref(z), acb_imagref(x), y, prec);
}

ACB_INLINE void
acb_addmul_arb(acb_t z, const acb_t x, const arb_t y, long prec)
{
    arb_addmul(acb_realref(z), acb_realref(x), y, prec);
    arb_addmul(acb_imagref(z), acb_imagref(x), y, prec);
}

ACB_INLINE void
acb_submul_arb(acb_t z, const acb_t x, const arb_t y, long prec)
{
    arb_submul(acb_realref(z), acb_realref(x), y, prec);
    arb_submul(acb_imagref(z), acb_imagref(x), y, prec);
}

void acb_inv(acb_t z, const acb_t x, long prec);

void acb_div(acb_t z, const acb_t x, const acb_t y, long prec);

ACB_INLINE void
acb_div_ui(acb_t z, const acb_t x, ulong c, long prec)
{
    arb_div_ui(acb_realref(z), acb_realref(x), c, prec);
    arb_div_ui(acb_imagref(z), acb_imagref(x), c, prec);
}

ACB_INLINE void
acb_div_si(acb_t z, const acb_t x, long c, long prec)
{
    arb_div_si(acb_realref(z), acb_realref(x), c, prec);
    arb_div_si(acb_imagref(z), acb_imagref(x), c, prec);
}

ACB_INLINE void
acb_div_arb(acb_t z, const acb_t x, const arb_t c, long prec)
{
    arb_div(acb_realref(z), acb_realref(x), c, prec);
    arb_div(acb_imagref(z), acb_imagref(x), c, prec);
}

ACB_INLINE void
acb_div_fmpz(acb_t z, const acb_t x, const fmpz_t c, long prec)
{
    arb_div_fmpz(acb_realref(z), acb_realref(x), c, prec);
    arb_div_fmpz(acb_imagref(z), acb_imagref(x), c, prec);
}

void acb_cube(acb_t y, const acb_t x, long prec);
void acb_pow_fmpz(acb_t y, const acb_t b, const fmpz_t e, long prec);
void acb_pow_ui(acb_t y, const acb_t b, ulong e, long prec);
void acb_pow_si(acb_t y, const acb_t b, long e, long prec);

ACB_INLINE void
acb_const_pi(acb_t x, long prec)
{
    arb_const_pi(acb_realref(x), prec);
    arb_zero(acb_imagref(x));
}

void acb_log(acb_t r, const acb_t z, long prec);
void acb_log1p(acb_t r, const acb_t z, long prec);
void acb_atan(acb_t r, const acb_t z, long prec);

void acb_exp(acb_t r, const acb_t z, long prec);
void acb_exp_pi_i(acb_t r, const acb_t z, long prec);

void acb_sin(acb_t r, const acb_t z, long prec);
void acb_cos(acb_t r, const acb_t z, long prec);
void acb_sin_cos(acb_t s, acb_t c, const acb_t z, long prec);
void acb_tan(acb_t r, const acb_t z, long prec);
void acb_cot(acb_t r, const acb_t z, long prec);

void acb_sin_pi(acb_t r, const acb_t z, long prec);
void acb_cos_pi(acb_t r, const acb_t z, long prec);
void acb_sin_cos_pi(acb_t s, acb_t c, const acb_t z, long prec);

void acb_tan_pi(acb_t r, const acb_t z, long prec);
void acb_cot_pi(acb_t r, const acb_t z, long prec);

void acb_pow_arb(acb_t z, const acb_t x, const arb_t y, long prec);
void acb_pow(acb_t r, const acb_t x, const acb_t y, long prec);

void acb_sqrt(acb_t y, const acb_t x, long prec);
void acb_rsqrt(acb_t y, const acb_t x, long prec);

void acb_rising_ui_bs(acb_t y, const acb_t x, ulong n, long prec);
void acb_rising_ui_rs(acb_t y, const acb_t x, ulong n, ulong m, long prec);
void acb_rising_ui_rec(acb_t y, const acb_t x, ulong n, long prec);
void acb_rising_ui(acb_t z, const acb_t x, ulong n, long prec);

void acb_rising2_ui_bs(acb_t u, acb_t v, const acb_t x, ulong n, long prec);
void acb_rising2_ui_rs(acb_t u, acb_t v, const acb_t x, ulong n, ulong m, long prec);
void acb_rising2_ui(acb_t u, acb_t v, const acb_t x, ulong n, long prec);

void acb_rising_ui_get_mag(mag_t bound, const acb_t s, ulong n);

void acb_gamma(acb_t y, const acb_t x, long prec);
void acb_rgamma(acb_t y, const acb_t x, long prec);
void acb_lgamma(acb_t y, const acb_t x, long prec);
void acb_digamma(acb_t y, const acb_t x, long prec);
void acb_zeta(acb_t z, const acb_t s, long prec);
void acb_hurwitz_zeta(acb_t z, const acb_t s, const acb_t a, long prec);

void acb_polylog(acb_t w, const acb_t s, const acb_t z, long prec);
void acb_polylog_si(acb_t w, long s, const acb_t z, long prec);

void acb_agm1(acb_t m, const acb_t z, long prec);
void acb_agm1_cpx(acb_ptr m, const acb_t z, long len, long prec);

/*
TBD

void acb_invroot_newton(acb_t r, const acb_t a, ulong m, const acb_t r0, long startprec, long prec);
void acb_root_exp(acb_t r, const acb_t a, long m, long index, long prec);
void acb_root_newton(acb_t r, const acb_t a, long m, long index, long prec);
void acb_root(acb_t r, const acb_t a, long m, long index, long prec);
*/

/* TODO: document */
ACB_INLINE int
acb_is_finite(const acb_t x)
{
    return arb_is_finite(acb_realref(x)) && arb_is_finite(acb_imagref(x));
}

/* TODO: document */
ACB_INLINE void
acb_indeterminate(acb_t x)
{
    arb_indeterminate(acb_realref(x));
    arb_indeterminate(acb_imagref(x));
}

ACB_INLINE void
_acb_vec_zero(acb_ptr A, long n)
{
    long i;
    for (i = 0; i < n; i++)
        acb_zero(A + i);
}

ACB_INLINE int
_acb_vec_is_zero(acb_srcptr vec, long len)
{
    long i;
    for (i = 0; i < len; i++)
        if (!acb_is_zero(vec + i))
            return 0;
    return 1;
}

ACB_INLINE void
_acb_vec_set(acb_ptr res, acb_srcptr vec, long len)
{
    long i;
    for (i = 0; i < len; i++)
        acb_set(res + i, vec + i);
}

ACB_INLINE void
_acb_vec_set_round(acb_ptr res, acb_srcptr vec, long len, long prec)
{
    long i;
    for (i = 0; i < len; i++)
        acb_set_round(res + i, vec + i, prec);
}

ACB_INLINE void
_acb_vec_neg(acb_ptr res, acb_srcptr vec, long len)
{
    long i;
    for (i = 0; i < len; i++)
        acb_neg(res + i, vec + i);
}

ACB_INLINE void
_acb_vec_add(acb_ptr res, acb_srcptr vec1, acb_srcptr vec2, long len, long prec)
{
    long i;
    for (i = 0; i < len; i++)
        acb_add(res + i, vec1 + i, vec2 + i, prec);
}

ACB_INLINE void
_acb_vec_sub(acb_ptr res, acb_srcptr vec1, acb_srcptr vec2, long len, long prec)
{
    long i;
    for (i = 0; i < len; i++)
        acb_sub(res + i, vec1 + i, vec2 + i, prec);
}

ACB_INLINE void
_acb_vec_scalar_submul(acb_ptr res, acb_srcptr vec, long len, const acb_t c, long prec)
{
    long i;
    for (i = 0; i < len; i++)
        acb_submul(res + i, vec + i, c, prec);
}

ACB_INLINE void
_acb_vec_scalar_addmul(acb_ptr res, acb_srcptr vec, long len, const acb_t c, long prec)
{
    long i;
    for (i = 0; i < len; i++)
        acb_addmul(res + i, vec + i, c, prec);
}

ACB_INLINE void
_acb_vec_scalar_mul(acb_ptr res, acb_srcptr vec, long len, const acb_t c, long prec)
{
    long i;
    for (i = 0; i < len; i++)
        acb_mul(res + i, vec + i, c, prec);
}

ACB_INLINE void
_acb_vec_scalar_mul_ui(acb_ptr res, acb_srcptr vec, long len, ulong c, long prec)
{
    long i;
    for (i = 0; i < len; i++)
        acb_mul_ui(res + i, vec + i, c, prec);
}

ACB_INLINE void
_acb_vec_scalar_mul_2exp_si(acb_ptr res, acb_srcptr vec, long len, long c)
{
    long i;
    for (i = 0; i < len; i++)
        acb_mul_2exp_si(res + i, vec + i, c);
}

ACB_INLINE void
_acb_vec_scalar_div_ui(acb_ptr res, acb_srcptr vec, long len, ulong c, long prec)
{
    long i;
    for (i = 0; i < len; i++)
        acb_div_ui(res + i, vec + i, c, prec);
}

ACB_INLINE void
_acb_vec_scalar_div(acb_ptr res, acb_srcptr vec, long len, const acb_t c, long prec)
{
    long i;
    for (i = 0; i < len; i++)
        acb_div(res + i, vec + i, c, prec);
}

ACB_INLINE void
_acb_vec_scalar_mul_arb(acb_ptr res, acb_srcptr vec, long len, const arb_t c, long prec)
{
    long i;
    for (i = 0; i < len; i++)
        acb_mul_arb(res + i, vec + i, c, prec);
}

ACB_INLINE void
_acb_vec_scalar_div_arb(acb_ptr res, acb_srcptr vec, long len, const arb_t c, long prec)
{
    long i;
    for (i = 0; i < len; i++)
    {
        arb_div(acb_realref(res + i), acb_realref(vec + i), c, prec);
        arb_div(acb_imagref(res + i), acb_imagref(vec + i), c, prec);
    }
}

ACB_INLINE void
_acb_vec_scalar_mul_fmpz(acb_ptr res, acb_srcptr vec, long len, const fmpz_t c, long prec)
{
    long i;
    for (i = 0; i < len; i++)
        acb_mul_fmpz(res + i, vec + i, c, prec);
}

ACB_INLINE void
_acb_vec_scalar_div_fmpz(acb_ptr res, acb_srcptr vec, long len, const fmpz_t c, long prec)
{
    long i;
    for (i = 0; i < len; i++)
        acb_div_fmpz(res + i, vec + i, c, prec);
}

ACB_INLINE void
acb_print(const acb_t x)
{
    printf("(");
    arb_print(acb_realref(x));
    printf(", ");
    arb_print(acb_imagref(x));
    printf(")");
}

void acb_printd(const acb_t z, long digits);

void acb_randtest(acb_t z, flint_rand_t state, long prec, long mag_bits);

void acb_randtest_special(acb_t z, flint_rand_t state, long prec, long mag_bits);

void acb_randtest_precise(acb_t z, flint_rand_t state, long prec, long mag_bits);

long acb_rel_error_bits(const acb_t x);

ARB_INLINE long
acb_rel_accuracy_bits(const acb_t x)
{
    return -acb_rel_error_bits(x);
}


ACB_INLINE long
acb_bits(const acb_t x)
{
    long b1, b2;
    b1 = arb_bits(acb_realref(x));
    b2 = arb_bits(acb_imagref(x));
    return FLINT_MAX(b1, b2);
}

ACB_INLINE int
acb_is_real(const acb_t x)
{
    return arb_is_zero(acb_imagref(x));
}

ACB_INLINE int
_acb_vec_is_real(acb_srcptr v, long len)
{
    long i;

    for (i = 0; i < len; i++)
    {
        if (!acb_is_real(v + i))
            return 0;
    }

    return 1;
}

ACB_INLINE long
_acb_vec_bits(acb_srcptr vec, long len)
{
    return _arb_vec_bits((arb_srcptr) vec, 2 * len);
}

ACB_INLINE void
_acb_vec_set_powers(acb_ptr xs, const acb_t x, long len, long prec)
{
    long i;

    for (i = 0; i < len; i++)
    {
        if (i == 0)
            acb_one(xs + i);
        else if (i == 1)
            acb_set_round(xs + i, x, prec);
        else if (i % 2 == 0)
            acb_mul(xs + i, xs + i / 2, xs + i / 2, prec);
        else
            acb_mul(xs + i, xs + i - 1, x, prec);
    }
}

ACB_INLINE void
_acb_vec_add_error_arf_vec(acb_ptr res, arf_srcptr err, long len)
{
    long i;
    for (i = 0; i < len; i++)
        acb_add_error_arf(res + i, err + i);
}

ACB_INLINE void
_acb_vec_add_error_mag_vec(acb_ptr res, mag_srcptr err, long len)
{
    long i;
    for (i = 0; i < len; i++)
    {
        mag_add(arb_radref(acb_realref(res + i)),
            arb_radref(acb_realref(res + i)), err + i);
        mag_add(arb_radref(acb_imagref(res + i)),
            arb_radref(acb_imagref(res + i)), err + i);
    }
}

ACB_INLINE void
_acb_vec_indeterminate(acb_ptr vec, long len)
{
    _arb_vec_indeterminate((arb_ptr) vec, 2 * len);
}

ACB_INLINE void
_acb_vec_trim(acb_ptr res, acb_srcptr vec, long len)
{
    long i;
    for (i = 0; i < len; i++)
        acb_trim(res + i, vec + i);
}

ARB_INLINE int
_acb_vec_get_unique_fmpz_vec(fmpz * res,  acb_srcptr vec, long len)
{
    long i;
    for (i = 0; i < len; i++)
        if (!acb_get_unique_fmpz(res + i, vec + i))
            return 0;
    return 1;
}

/* sort complex numbers in a nice-to-display order */
void _acb_vec_sort_pretty(acb_ptr vec, long len);

#ifdef __cplusplus
}
#endif

#endif
