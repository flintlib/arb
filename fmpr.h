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
#include <mpir.h>
#include <mpfr.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpq.h"

#define fmpr_rnd_t int
#define FMPR_RND_FLOOR 0
#define FMPR_RND_CEIL 1
#define FMPR_RND_UP 2
#define FMPR_RND_DOWN 3
#define FMPR_RND_NEAR 4

typedef struct
{
    fmpz man;
    fmpz exp;
}
fmpr_struct;

typedef fmpr_struct fmpr_t[1];

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

/* The return value encodes an error bound, or FMPR_RESULT_EXACT if exact */

#define FMPR_RESULT_EXACT LONG_MAX
#define FMPR_PREC_EXACT LONG_MAX

/* TODO: version for just the val2 reduction! */
long _fmpr_normalise(fmpz_t man, fmpz_t exp, long prec, fmpr_rnd_t rnd);

static __inline__ void
fmpr_set(fmpr_t y, const fmpr_t x)
{
    if (y != x)
    {
        fmpz_set(fmpr_manref(y), fmpr_manref(x));
        fmpz_set(fmpr_expref(y), fmpr_expref(x));
    }
}

static __inline__ long
fmpr_set_round(fmpr_t y, const fmpr_t x, long prec, fmpr_rnd_t rnd)
{
    fmpr_set(y, x);
    if (fmpr_is_special(y))
        return FMPR_RESULT_EXACT;
    else
        return _fmpr_normalise(fmpr_manref(y), fmpr_expref(y), prec, rnd);
}


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




void fmpr_randtest(fmpr_t x, flint_rand_t state, long bits, long exp_bits);

void fmpr_randtest_not_zero(fmpr_t x, flint_rand_t state, long bits, long exp_bits);

void fmpr_randtest_special(fmpr_t x, flint_rand_t state, long bits, long exp_bits);

static __inline__ int
fmpr_get_mpfr(mpfr_t x, const fmpr_t y, mpfr_rnd_t rnd)
{
    int r;

    if (fmpr_is_special(y))
    {
        if (fmpr_is_zero(y)) mpfr_set_zero(x, 0);
        else if (fmpr_is_pos_inf(y)) mpfr_set_inf(x, 1);
        else if (fmpr_is_neg_inf(y)) mpfr_set_inf(x, -1);
        else mpfr_set_nan(x);
        r = 0;
    }
    else
    {
        if (!COEFF_IS_MPZ(*fmpr_manref(y)))
            r = mpfr_set_si_2exp(x, *fmpr_manref(y), *fmpr_expref(y), rnd);
        else
            r = mpfr_set_z_2exp(x, COEFF_TO_PTR(*fmpr_manref(y)), *fmpr_expref(y), rnd);
    }

    return r;
}

static __inline__ void
fmpr_set_mpfr(fmpr_t x, const mpfr_t y)
{
    if (!mpfr_regular_p(y))
    {
        if (mpfr_zero_p(y)) fmpr_zero(x);
        else if (mpfr_inf_p(y) && mpfr_sgn(y) > 0) fmpr_pos_inf(x);
        else if (mpfr_inf_p(y) && mpfr_sgn(y) < 0) fmpr_neg_inf(x);
        else fmpr_nan(x);
    }
    else
    {
        __mpz_struct * z = _fmpz_promote(fmpr_manref(x));
        fmpz_set_si(fmpr_expref(x), mpfr_get_z_2exp(z, y));
        _fmpz_demote_val(fmpr_manref(x));
        _fmpr_normalise(fmpr_manref(x), fmpr_expref(x), y->_mpfr_prec, FMPR_RND_DOWN);
    }
}

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

/* Some fmpz helper functions */

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

static long _fmpz_sub_small_large(const fmpz_t x, const fmpz_t y)
{
    fmpz_t t;
    fmpz_init(t);
    fmpz_sub(t, x, y);
    if (!COEFF_IS_MPZ(*t))
    {
        /* no need to free t */
        return *t;
    }
    else
    {
        int sign = fmpz_sgn(t);
        fmpz_clear(t);
        return (sign > 0) ? LONG_MAX : -LONG_MAX;
    }
}

static __inline__ long _fmpz_sub_small(const fmpz_t x, const fmpz_t y)
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


/* Arithmetic */


long _fmpr_add_eps(fmpr_t z, const fmpr_t x, int sign, long prec, fmpr_rnd_t rnd);

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

long fmpr_div(fmpr_t z, const fmpr_t x, const fmpr_t y, long prec, fmpr_rnd_t rnd);
long fmpr_div_ui(fmpr_t z, const fmpr_t x, ulong y, long prec, fmpr_rnd_t rnd);
long fmpr_ui_div(fmpr_t z, ulong x, const fmpr_t y, long prec, fmpr_rnd_t rnd);
long fmpr_div_si(fmpr_t z, const fmpr_t x, long y, long prec, fmpr_rnd_t rnd);
long fmpr_si_div(fmpr_t z, long x, const fmpr_t y, long prec, fmpr_rnd_t rnd);
long fmpr_div_fmpz(fmpr_t z, const fmpr_t x, const fmpz_t y, long prec, fmpr_rnd_t rnd);
long fmpr_fmpz_div(fmpr_t z, const fmpz_t x, const fmpr_t y, long prec, fmpr_rnd_t rnd);
long fmpr_fmpz_div_fmpz(fmpr_t z, const fmpz_t x, const fmpz_t y, long prec, fmpr_rnd_t rnd);

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

long fmpr_log(fmpr_t y, const fmpr_t x, long prec, fmpr_rnd_t rnd);

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
        if (e > 0)
            fmpz_add_ui(fmpr_expref(y), fmpr_expref(x), e);
        else
            fmpz_sub_ui(fmpr_expref(y), fmpr_expref(x), -e);
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

#endif

