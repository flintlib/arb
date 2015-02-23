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

#ifndef FMPR_H
#define FMPR_H

#include <stdio.h>
#include <limits.h>
#include <gmp.h>
#include <mpfr.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpq.h"
#include "fmpz_extras.h"

#include "config.h"
#ifdef HAVE_TLS
#if HAVE_TLS
#define TLS_PREFIX __thread
#else
#define TLS_PREFIX
#endif
#else
#define TLS_PREFIX
#endif

#define fmpr_rnd_t int
#define FMPR_RND_DOWN 0
#define FMPR_RND_UP 1
#define FMPR_RND_FLOOR 2
#define FMPR_RND_CEIL 3
#define FMPR_RND_NEAR 4

#ifdef __cplusplus
extern "C" {
#endif

static __inline__ int
rounds_up(fmpr_rnd_t rnd, int negative)
{
    if (rnd == FMPR_RND_DOWN) return 0;
    if (rnd == FMPR_RND_UP) return 1;
    if (rnd == FMPR_RND_FLOOR) return negative;
    return !negative;
}

static __inline__ mpfr_rnd_t rnd_to_mpfr(fmpr_rnd_t rnd)
{
    if (rnd == FMPR_RND_DOWN) return MPFR_RNDZ;
    else if (rnd == FMPR_RND_UP) return MPFR_RNDA;
    else if (rnd == FMPR_RND_FLOOR) return MPFR_RNDD;
    else if (rnd == FMPR_RND_CEIL) return MPFR_RNDU;
    else return MPFR_RNDN;
}

typedef struct
{
    fmpz man;
    fmpz exp;
}
fmpr_struct;

typedef fmpr_struct fmpr_t[1];
typedef fmpr_struct * fmpr_ptr;
typedef const fmpr_struct * fmpr_srcptr;

#define fmpr_manref(x) (&(x)->man)
#define fmpr_expref(x) (&(x)->exp)


static __inline__ void
fmpr_init(fmpr_t x)
{
    fmpz_init(fmpr_manref(x));
    fmpz_init(fmpr_expref(x));
}

static __inline__ void
fmpr_clear(fmpr_t x)
{
    fmpz_clear(fmpr_manref(x));
    fmpz_clear(fmpr_expref(x));
}


/* The return value encodes an error bound, or FMPR_RESULT_EXACT if exact */
#define FMPR_RESULT_EXACT LONG_MAX

/* Allow 'infinite' precision for operations where a result can be computed exactly */
#define FMPR_PREC_EXACT LONG_MAX

#define FMPR_PREC_ADD(prec,extra) ((prec) == FMPR_PREC_EXACT ? FMPR_PREC_EXACT : (prec) + (extra))

/* ------------------------------------------------------------------------ */
/* Special values */

/*
  Special values (0, +inf, -inf, NaN) are encoded using a zero mantissa.
  Zero is encoded with a zero exponent.
*/

#define FMPR_EXP_POS_INF 1L
#define FMPR_EXP_NEG_INF 2L
#define FMPR_EXP_NAN 3L

static __inline__ void fmpr_zero(fmpr_t x) {
    fmpz_zero(fmpr_manref(x));
    fmpz_zero(fmpr_expref(x));
}

static __inline__ int fmpr_is_zero(const fmpr_t x) {
    return fmpz_is_zero(fmpr_manref(x)) && fmpz_is_zero(fmpr_expref(x));
}

static __inline__ int fmpr_is_special(const fmpr_t x) {
    return fmpz_is_zero(fmpr_manref(x));
}

static __inline__ int fmpr_is_normal(const fmpr_t x) {
    return !fmpz_is_zero(fmpr_manref(x));
}

static __inline__ int fmpr_is_inf(const fmpr_t x) {
    return fmpr_is_special(x) && (*fmpr_expref(x) == FMPR_EXP_POS_INF ||
                                  *fmpr_expref(x) == FMPR_EXP_NEG_INF);
}

static __inline__ int fmpr_is_pos_inf(const fmpr_t x) {
    return fmpr_is_special(x) && (*fmpr_expref(x) == FMPR_EXP_POS_INF);
}

static __inline__ int fmpr_is_neg_inf(const fmpr_t x) {
    return fmpr_is_special(x) && (*fmpr_expref(x) == FMPR_EXP_NEG_INF);
}

static __inline__ int fmpr_is_nan(const fmpr_t x) {
    return fmpr_is_special(x) && (*fmpr_expref(x) == FMPR_EXP_NAN);
}

static __inline__ int fmpr_is_finite(const fmpr_t x) {
    return fmpr_is_zero(x) || !fmpr_is_special(x);
}

static __inline__ void fmpr_pos_inf(fmpr_t x) {
    fmpz_zero(fmpr_manref(x));
    fmpz_set_si(fmpr_expref(x), FMPR_EXP_POS_INF);
}

static __inline__ void fmpr_neg_inf(fmpr_t x) {
    fmpz_zero(fmpr_manref(x));
    fmpz_set_si(fmpr_expref(x), FMPR_EXP_NEG_INF);
}

static __inline__ void fmpr_nan(fmpr_t x) {
    fmpz_zero(fmpr_manref(x));
    fmpz_set_si(fmpr_expref(x), FMPR_EXP_NAN);
}



static __inline__ int
fmpr_is_one(const fmpr_t x)
{
    return fmpz_is_one(fmpr_manref(x)) && fmpz_is_zero(fmpr_expref(x));
}

static __inline__ void
fmpr_one(fmpr_t x)
{
    fmpz_one(fmpr_manref(x));
    fmpz_zero(fmpr_expref(x));
}

/* ------------------------------------------------------------------------ */

long _fmpr_normalise_naive(fmpz_t man, fmpz_t exp, long prec, fmpr_rnd_t rnd);

static __inline__ void
fmpr_set(fmpr_t y, const fmpr_t x)
{
    if (y != x)
    {
        fmpz_set(fmpr_manref(y), fmpr_manref(x));
        fmpz_set(fmpr_expref(y), fmpr_expref(x));
    }
}

static __inline__ void
fmpr_swap(fmpr_t x, fmpr_t y)
{
    fmpz_swap(fmpr_manref(x), fmpr_manref(y));
    fmpz_swap(fmpr_expref(x), fmpr_expref(y));
}

long _fmpr_set_round(fmpz_t rman, fmpz_t rexp,
    const fmpz_t man, const fmpz_t exp, long prec, fmpr_rnd_t rnd);

static __inline__ long
_fmpr_normalise(fmpz_t man, fmpz_t exp, long prec, fmpr_rnd_t rnd)
{
    if (fmpz_is_zero(man))
    {
        fmpz_zero(man);
        fmpz_zero(exp);
        return FMPR_RESULT_EXACT;
    }
    else
    {
        return _fmpr_set_round(man, exp, man, exp, prec, rnd);
    }
}

static __inline__ long
fmpr_set_round_naive(fmpr_t y, const fmpr_t x, long prec, fmpr_rnd_t rnd)
{
    fmpr_set(y, x);
    if (fmpr_is_special(y))
        return FMPR_RESULT_EXACT;
    else
        return _fmpr_normalise(fmpr_manref(y), fmpr_expref(y), prec, rnd);
}

static __inline__ long
fmpr_set_round(fmpr_t y, const fmpr_t x, long prec, fmpr_rnd_t rnd)
{
    if (fmpr_is_special(x))
    {
        fmpr_set(y, x);
        return FMPR_RESULT_EXACT;
    }
    else
    {
        return _fmpr_set_round(fmpr_manref(y), fmpr_expref(y),
            fmpr_manref(x), fmpr_expref(x), prec, rnd);
    }
}

static __inline__ long
fmpr_set_round_fmpz_2exp(fmpr_t y, const fmpz_t x, const fmpz_t exp, long prec, fmpr_rnd_t rnd)
{
    if (fmpz_is_zero(x))
    {
        fmpr_zero(y);
        return FMPR_RESULT_EXACT;
    }
    else
    {
        return _fmpr_set_round(fmpr_manref(y), fmpr_expref(y), x, exp, prec, rnd);
    }
}

static __inline__ long
fmpr_set_round_fmpz(fmpr_t y, const fmpz_t x, long prec, fmpr_rnd_t rnd)
{
    if (fmpz_is_zero(x))
    {
        fmpr_zero(y);
        return FMPR_RESULT_EXACT;
    }
    else
    {
        long ret;
        fmpz_t exp;
        fmpz_init(exp);
        ret = _fmpr_set_round(fmpr_manref(y), fmpr_expref(y), x, exp, prec, rnd);
        fmpz_clear(exp);
        return ret;
    }
}

long _fmpr_set_round_mpn(long * shift, fmpz_t man, mp_srcptr x, mp_size_t xn, int negative, long prec, fmpr_rnd_t rnd);

long fmpr_set_round_ui_2exp_fmpz(fmpr_t z,
        mp_limb_t lo, const fmpz_t exp, int negative,
        long prec, fmpr_rnd_t rnd);

long fmpr_set_round_uiui_2exp_fmpz(fmpr_t z,
    mp_limb_t hi, mp_limb_t lo, const fmpz_t exp, int negative,
    long prec, fmpr_rnd_t rnd);

void fmpr_ulp(fmpr_t u, const fmpr_t x, long prec);

int fmpr_check_ulp(const fmpr_t result, long r, long prec);

static __inline__ int
fmpr_equal(const fmpr_t x, const fmpr_t y)
{
    return fmpz_equal(fmpr_expref(x), fmpr_expref(y)) &&
        fmpz_equal(fmpr_manref(x), fmpr_manref(y));
}

static __inline__ int
fmpr_sgn(const fmpr_t x)
{
    if (fmpr_is_special(x))
    {
        if (fmpr_is_zero(x))
            return 0;
        if (fmpr_is_pos_inf(x))
            return 1;
        if (fmpr_is_neg_inf(x))
            return -1;
        return 0;
    }

    return fmpz_sgn(fmpr_manref(x));
}

int fmpr_cmp(const fmpr_t x, const fmpr_t y);

int fmpr_cmpabs(const fmpr_t x, const fmpr_t y);

int fmpr_cmpabs_ui(const fmpr_t x, ulong y);

void fmpr_randtest(fmpr_t x, flint_rand_t state, long bits, long exp_bits);

void fmpr_randtest_not_zero(fmpr_t x, flint_rand_t state, long bits, long exp_bits);

void fmpr_randtest_special(fmpr_t x, flint_rand_t state, long bits, long exp_bits);

int fmpr_get_mpfr(mpfr_t x, const fmpr_t y, mpfr_rnd_t rnd);

void fmpr_set_mpfr(fmpr_t x, const mpfr_t y);

double fmpr_get_d(const fmpr_t x, fmpr_rnd_t rnd);

void fmpr_set_d(fmpr_t x, double v);

static __inline__ void fmpr_set_ui(fmpr_t x, ulong c)
{
    if (c == 0)
    {
        fmpr_zero(x);
    }
    else
    {
        int b;
        count_trailing_zeros(b, c);
        fmpz_set_ui(fmpr_manref(x), c >> b);
        fmpz_set_ui(fmpr_expref(x), b);
    }
}

static __inline__ void fmpr_set_si(fmpr_t x, long c)
{
    if (c == 0)
    {
        fmpr_zero(x);
    }
    else
    {
        int b;
        count_trailing_zeros(b, c);
        fmpz_set_si(fmpr_manref(x), c >> b);
        fmpz_set_ui(fmpr_expref(x), b);
    }
}

static __inline__ void
fmpr_set_fmpz(fmpr_t x, const fmpz_t c)
{
    if (fmpz_is_zero(c))
    {
        fmpr_zero(x);
    }
    else
    {
        long v = fmpz_val2(c);

        fmpz_tdiv_q_2exp(fmpr_manref(x), c, v);
        fmpz_set_ui(fmpr_expref(x), v);
    }
}

/* Arithmetic */


long _fmpr_add_eps(fmpr_t z, const fmpr_t x, int sign, long prec, fmpr_rnd_t rnd);

long _fmpr_add_mpn(fmpr_t z,
        mp_srcptr xman, mp_size_t xn, int xsign, const fmpz_t xexp,
        mp_srcptr yman, mp_size_t yn, int ysign, const fmpz_t yexp,
        long shift, long prec, fmpr_rnd_t rnd);

static __inline__ long
_fmpr_add_1x1(fmpr_t z,
        mp_limb_t x, int xsign, const fmpz_t xexp,
        mp_limb_t y, int ysign, const fmpz_t yexp,
        long shift, long prec, long rnd)
{
    mp_limb_t hi, lo, t, u;
    int sign = ysign;

    t = x;
    u = y;

    (void) yexp; /* unused */

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

long fmpr_add_naive(fmpr_t z, const fmpr_t x, const fmpr_t y, long prec, fmpr_rnd_t rnd);

long _fmpr_mul_mpn(fmpr_t z,
    mp_srcptr xman, mp_size_t xn, const fmpz_t xexp,
    mp_srcptr yman, mp_size_t yn, const fmpz_t yexp,
    int negative, long prec, fmpr_rnd_t rnd);

long _fmpr_mul_1x1(fmpr_t z,
  mp_limb_t u, const fmpz_t xexp, mp_limb_t v, const fmpz_t yexp,
  int negative, long prec, fmpr_rnd_t rnd);

long fmpr_mul_naive(fmpr_t z, const fmpr_t x, const fmpr_t y, long prec, fmpr_rnd_t rnd);
long fmpr_mul(fmpr_t z, const fmpr_t x, const fmpr_t y, long prec, fmpr_rnd_t rnd);
long fmpr_mul_ui(fmpr_t z, const fmpr_t x, ulong y, long prec, fmpr_rnd_t rnd);
long fmpr_mul_si(fmpr_t z, const fmpr_t x, long y, long prec, fmpr_rnd_t rnd);
long fmpr_mul_fmpz(fmpr_t z, const fmpr_t x, const fmpz_t y, long prec, fmpr_rnd_t rnd);

long fmpr_add(fmpr_t z, const fmpr_t x, const fmpr_t y, long prec, fmpr_rnd_t rnd);
long fmpr_add_ui(fmpr_t z, const fmpr_t x, ulong y, long prec, fmpr_rnd_t rnd);
long fmpr_add_si(fmpr_t z, const fmpr_t x, long y, long prec, fmpr_rnd_t rnd);
long fmpr_add_fmpz(fmpr_t z, const fmpr_t x, const fmpz_t y, long prec, fmpr_rnd_t rnd);

long fmpr_sub(fmpr_t z, const fmpr_t x, const fmpr_t y, long prec, fmpr_rnd_t rnd);
long fmpr_sub_ui(fmpr_t z, const fmpr_t x, ulong y, long prec, fmpr_rnd_t rnd);
long fmpr_sub_si(fmpr_t z, const fmpr_t x, long y, long prec, fmpr_rnd_t rnd);
long fmpr_sub_fmpz(fmpr_t z, const fmpr_t x, const fmpz_t y, long prec, fmpr_rnd_t rnd);

long fmpr_sum(fmpr_t s, const fmpr_struct * terms, long len, long prec, fmpr_rnd_t rnd);

long fmpr_div(fmpr_t z, const fmpr_t x, const fmpr_t y, long prec, fmpr_rnd_t rnd);
long fmpr_div_ui(fmpr_t z, const fmpr_t x, ulong y, long prec, fmpr_rnd_t rnd);
long fmpr_ui_div(fmpr_t z, ulong x, const fmpr_t y, long prec, fmpr_rnd_t rnd);
long fmpr_div_si(fmpr_t z, const fmpr_t x, long y, long prec, fmpr_rnd_t rnd);
long fmpr_si_div(fmpr_t z, long x, const fmpr_t y, long prec, fmpr_rnd_t rnd);
long fmpr_div_fmpz(fmpr_t z, const fmpr_t x, const fmpz_t y, long prec, fmpr_rnd_t rnd);
long fmpr_fmpz_div(fmpr_t z, const fmpz_t x, const fmpr_t y, long prec, fmpr_rnd_t rnd);
long fmpr_fmpz_div_fmpz(fmpr_t z, const fmpz_t x, const fmpz_t y, long prec, fmpr_rnd_t rnd);

void fmpr_divappr_abs_ubound(fmpr_t z, const fmpr_t x, const fmpr_t y, long prec);

long fmpr_addmul(fmpr_t z, const fmpr_t x, const fmpr_t y, long prec, fmpr_rnd_t rnd);
long fmpr_addmul_ui(fmpr_t z, const fmpr_t x, ulong y, long prec, fmpr_rnd_t rnd);
long fmpr_addmul_si(fmpr_t z, const fmpr_t x, long y, long prec, fmpr_rnd_t rnd);
long fmpr_addmul_fmpz(fmpr_t z, const fmpr_t x, const fmpz_t y, long prec, fmpr_rnd_t rnd);

long fmpr_submul(fmpr_t z, const fmpr_t x, const fmpr_t y, long prec, fmpr_rnd_t rnd);
long fmpr_submul_ui(fmpr_t z, const fmpr_t x, ulong y, long prec, fmpr_rnd_t rnd);
long fmpr_submul_si(fmpr_t z, const fmpr_t x, long y, long prec, fmpr_rnd_t rnd);
long fmpr_submul_fmpz(fmpr_t z, const fmpr_t x, const fmpz_t y, long prec, fmpr_rnd_t rnd);

long fmpr_sqrt(fmpr_t y, const fmpr_t x, long prec, fmpr_rnd_t rnd);
long fmpr_sqrt_ui(fmpr_t z, ulong x, long prec, fmpr_rnd_t rnd);
long fmpr_sqrt_fmpz(fmpr_t z, const fmpz_t x, long prec, fmpr_rnd_t rnd);

long fmpr_rsqrt(fmpr_t y, const fmpr_t x, long prec, fmpr_rnd_t rnd);

long fmpr_root(fmpr_t y, const fmpr_t x, ulong k, long prec, fmpr_rnd_t rnd);

long fmpr_log(fmpr_t y, const fmpr_t x, long prec, fmpr_rnd_t rnd);
long fmpr_log1p(fmpr_t y, const fmpr_t x, long prec, fmpr_rnd_t rnd);

long fmpr_exp(fmpr_t y, const fmpr_t x, long prec, fmpr_rnd_t rnd);
long fmpr_expm1(fmpr_t y, const fmpr_t x, long prec, fmpr_rnd_t rnd);

static __inline__ void
fmpr_neg(fmpr_t y, const fmpr_t x)
{
    if (fmpr_is_special(x))
    {
        if (fmpr_is_pos_inf(x))
            fmpr_neg_inf(y);
        else if (fmpr_is_neg_inf(x))
            fmpr_pos_inf(y);
        else
            fmpr_set(y, x);
    }
    else
    {
        fmpz_neg(fmpr_manref(y), fmpr_manref(x));
        fmpz_set(fmpr_expref(y), fmpr_expref(x));
    }
}

static __inline__ long
fmpr_neg_round(fmpr_t y, const fmpr_t x, long prec, fmpr_rnd_t rnd)
{
    fmpr_neg(y, x);
    if (fmpr_is_special(y))
        return FMPR_RESULT_EXACT;
    else
        return _fmpr_normalise(fmpr_manref(y), fmpr_expref(y), prec, rnd);
}

static __inline__ void
fmpr_abs(fmpr_t y, const fmpr_t x)
{
    if (fmpr_sgn(x) < 0)
        fmpr_neg(y, x);
    else
        fmpr_set(y, x);
}

/* Convert return value to error bound */
static __inline__ void
fmpr_set_error_result(fmpr_t err, const fmpr_t result, long rret)
{
    if (rret == FMPR_RESULT_EXACT)
    {
        fmpr_zero(err);
    }
    /*
        XXX: normally, a special value should not result from inexact arithmetic
        else if (fmpr_is_special(result))
        {
            fmpr_pos_inf(err);
        }
    */
    else
    {
        /* TODO: inline this */
        fmpz_sub_ui(fmpr_expref(err), fmpr_expref(result), rret);
        fmpz_one(fmpr_manref(err));
    }
}

static __inline__ void
fmpr_add_error_result(fmpr_t err, const fmpr_t err_in,
    const fmpr_t result, long rret, long prec, fmpr_rnd_t rnd)
{
    if (rret == FMPR_RESULT_EXACT)
    {
        fmpr_set(err, err_in);
    }
    else
    {
        fmpr_t t;
        fmpr_init(t);
        fmpr_set_error_result(t, result, rret);
        fmpr_add(err, err_in, t, prec, rnd);
        fmpr_clear(t);
    }
}

void fmpr_print(const fmpr_t x);

void fmpr_printd(const fmpr_t x, long digits);


static __inline__ void
fmpr_mul_2exp_si(fmpr_t y, const fmpr_t x, long e)
{
    if (e == 0 || fmpr_is_special(x))
    {
        fmpr_set(y, x);
    }
    else
    {
        fmpz_set(fmpr_manref(y), fmpr_manref(x));
        fmpz_add_si_inline(fmpr_expref(y), fmpr_expref(x), e);
    }
}

static __inline__ void
fmpr_mul_2exp_fmpz(fmpr_t y, const fmpr_t x, const fmpz_t e)
{
    if (e == 0 || fmpr_is_special(x))
    {
        fmpr_set(y, x);
    }
    else
    {
        fmpz_set(fmpr_manref(y), fmpr_manref(x));
        fmpz_add(fmpr_expref(y), fmpr_expref(x), e);
    }
}

void fmpr_get_fmpq(fmpq_t y, const fmpr_t x);

long fmpr_set_fmpq(fmpr_t x, const fmpq_t y, long prec, fmpr_rnd_t rnd);

void fmpr_get_fmpz(fmpz_t z, const fmpr_t x, fmpr_rnd_t rnd);

long fmpr_get_si(const fmpr_t x, fmpr_rnd_t rnd);

void fmpr_set_fmpz_2exp(fmpr_t x, const fmpz_t man, const fmpz_t exp);

void fmpr_get_fmpz_2exp(fmpz_t man, fmpz_t exp, const fmpr_t x);

int fmpr_get_fmpz_fixed_fmpz(fmpz_t y, const fmpr_t x, const fmpz_t e);

int fmpr_get_fmpz_fixed_si(fmpz_t y, const fmpr_t x, long e);

static __inline__ void
fmpr_set_si_2exp_si(fmpr_t x, long man, long exp)
{
    fmpr_set_si(x, man);
    fmpr_mul_2exp_si(x, x, exp);
}

static __inline__ void
fmpr_set_ui_2exp_si(fmpr_t x, ulong man, long exp)
{
    fmpr_set_ui(x, man);
    fmpr_mul_2exp_si(x, x, exp);
}

int fmpr_cmp_2exp_si(const fmpr_t x, long e);

int fmpr_cmpabs_2exp_si(const fmpr_t x, long e);

static __inline__ void
fmpr_abs_bound_le_2exp_fmpz(fmpz_t b, const fmpr_t x)
{
    if (fmpz_is_pm1(fmpr_manref(x)))
        fmpz_set(b, fmpr_expref(x));
    else
        fmpz_add_ui(b, fmpr_expref(x), fmpz_bits(fmpr_manref(x)));
}

static __inline__ void
fmpr_abs_bound_lt_2exp_fmpz(fmpz_t b, const fmpr_t x)
{
    fmpz_add_ui(b, fmpr_expref(x), fmpz_bits(fmpr_manref(x)));
}

long fmpr_abs_bound_lt_2exp_si(const fmpr_t x);

static __inline__ void
fmpr_min(fmpr_t z, const fmpr_t a, const fmpr_t b)
{
    if (fmpr_cmp(a, b) <= 0)
        fmpr_set(z, a);
    else
        fmpr_set(z, b);
}

static __inline__ void
fmpr_max(fmpr_t z, const fmpr_t a, const fmpr_t b)
{
    if (fmpr_cmp(a, b) > 0)
        fmpr_set(z, a);
    else
        fmpr_set(z, b);
}

static __inline__ long
fmpr_bits(const fmpr_t x)
{
    if (fmpr_is_special(x))
        return 0;
    else
        return fmpz_bits(fmpr_manref(x));
}

static __inline__ int
fmpr_is_int(const fmpr_t x)
{
    if (fmpr_is_special(x))
        return fmpr_is_zero(x);
    else
        return fmpz_sgn(fmpr_expref(x)) >= 0;
}

static __inline__ int
fmpr_is_int_2exp_si(const fmpr_t x, long e)
{
    if (fmpr_is_special(x))
        return fmpr_is_zero(x);
    else
        return fmpz_cmp_si(fmpr_expref(x), e) >= 0;
}

#define CALL_MPFR_FUNC(r, func, y, x, prec, rnd) \
    do { \
        mpfr_t __t, __u; \
        mpfr_rnd_t __rnd; \
        __rnd = rnd_to_mpfr(rnd); \
        mpfr_init2(__t, 2 + fmpz_bits(fmpr_manref(x))); \
        mpfr_init2(__u, FLINT_MAX(2, prec)); \
        mpfr_set_emin(MPFR_EMIN_MIN); \
        mpfr_set_emax(MPFR_EMAX_MAX); \
        fmpr_get_mpfr(__t, x, MPFR_RNDD); \
        r = func(__u, __t, __rnd); \
        if (mpfr_overflow_p() || mpfr_underflow_p()) \
        { \
            printf("exception: mpfr overflow\n"); \
            abort(); \
        } \
        fmpr_set_mpfr(y, __u); \
        if (r == 0) \
            r = FMPR_RESULT_EXACT; \
        else \
            r = prec - fmpz_bits(fmpr_manref(y)); \
        mpfr_clear(__t); \
        mpfr_clear(__u); \
    } while (0);


#define CALL_MPFR_FUNC_K(r, func, y, x, k, prec, rnd) \
    do { \
        mpfr_t __t, __u; \
        mpfr_rnd_t __rnd; \
        __rnd = rnd_to_mpfr(rnd); \
        mpfr_init2(__t, 2 + fmpz_bits(fmpr_manref(x))); \
        mpfr_init2(__u, FLINT_MAX(2, prec)); \
        mpfr_set_emin(MPFR_EMIN_MIN); \
        mpfr_set_emax(MPFR_EMAX_MAX); \
        fmpr_get_mpfr(__t, x, MPFR_RNDD); \
        r = func(__u, __t, k, __rnd); \
        if (mpfr_overflow_p() || mpfr_underflow_p()) \
        { \
            printf("exception: mpfr overflow\n"); \
            abort(); \
        } \
        fmpr_set_mpfr(y, __u); \
        if (r == 0) \
            r = FMPR_RESULT_EXACT; \
        else \
            r = prec - fmpz_bits(fmpr_manref(y)); \
        mpfr_clear(__t); \
        mpfr_clear(__u); \
    } while (0);

#define CALL_MPFR_FUNC_2X1(r1, r2, func, y1, y2, x, prec, rnd) \
    do { \
        mpfr_t __t, __u, __v; \
        mpfr_rnd_t __rnd; \
        __rnd = rnd_to_mpfr(rnd); \
        mpfr_init2(__t, 2 + fmpz_bits(fmpr_manref(x))); \
        mpfr_init2(__u, FLINT_MAX(2, prec)); \
        mpfr_init2(__v, FLINT_MAX(2, prec)); \
        mpfr_set_emin(MPFR_EMIN_MIN); \
        mpfr_set_emax(MPFR_EMAX_MAX); \
        fmpr_get_mpfr(__t, x, MPFR_RNDD); \
        func(__u, __v, __t, __rnd); \
        if (mpfr_overflow_p() || mpfr_underflow_p()) \
        { \
            printf("exception: mpfr overflow\n"); \
            abort(); \
        } \
        fmpr_set_mpfr(y1, __u); \
        r1 = prec - fmpz_bits(fmpr_manref(y1)); \
        fmpr_set_mpfr(y2, __v); \
        r2 = prec - fmpz_bits(fmpr_manref(y2)); \
        mpfr_clear(__t); \
        mpfr_clear(__u); \
        mpfr_clear(__v); \
    } while (0);

void fmpr_pow_sloppy_fmpz(fmpr_t y, const fmpr_t b, const fmpz_t e,
    long prec, fmpr_rnd_t rnd);

void fmpr_pow_sloppy_ui(fmpr_t y, const fmpr_t b, ulong e,
    long prec, fmpr_rnd_t rnd);

void fmpr_pow_sloppy_si(fmpr_t y, const fmpr_t b, long e,
    long prec, fmpr_rnd_t rnd);

/* vector functions */

static __inline__ fmpr_ptr
_fmpr_vec_init(long n)
{
    long i;
    fmpr_ptr v = (fmpr_ptr) flint_malloc(sizeof(fmpr_struct) * n);

    for (i = 0; i < n; i++)
        fmpr_init(v + i);

    return v;
}

static __inline__ void
_fmpr_vec_clear(fmpr_ptr v, long n)
{
    long i;
    for (i = 0; i < n; i++)
        fmpr_clear(v + i);
    flint_free(v);
}


#ifdef __cplusplus
}
#endif

#endif

