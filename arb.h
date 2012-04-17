/*=============================================================================

    This file is part of ARB.

    ARB is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    ARB is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with ARB; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#ifndef ARB_H
#define ARB_H

#include <stdio.h>
#include <mpir.h>
#include <mpfr.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpq.h"
#include "ufloat.h"


typedef struct
{
    fmpz mid;
    fmpz exp;
    fmpz rad;
    long prec;
}
arb_struct;

#define arb_midref(x) (&(x)->mid)
#define arb_expref(x) (&(x)->exp)
#define arb_radref(x) (&(x)->rad)
#define arb_prec(x) ((x)->prec)

typedef arb_struct arb_t[1];

void arb_init(arb_t x, long prec);
void arb_clear(arb_t x);
void arb_set(arb_t x, const arb_t y);

void arb_randtest(arb_t x, flint_rand_t state, long exp_bits);
void arb_debug(const arb_t x);

void arb_set_fmpz(arb_t x, const fmpz_t c);
void arb_set_si(arb_t x, long c);
void arb_set_ui(arb_t x, ulong c);

int arb_contains_zero(const arb_t x);
int arb_contains_ui(const arb_t x, ulong z);
int _arb_contains_fmpq(const fmpz_t mid, const fmpz_t rad, long exp,
    const fmpz_t num, const fmpz_t den);
int arb_contains_fmpq(const arb_t x, const fmpq_t q);
int arb_contains_fmpz(const arb_t x, const fmpz_t z);
int arb_contains_mpfr(const arb_t x, const mpfr_t z);

void _arb_get_rand_fmpq(fmpz_t num, flint_rand_t state, const fmpz_t den,
    const fmpz_t mid, const fmpz_t rad, long exp);
void arb_get_rand_fmpq(fmpq_t q, flint_rand_t state, const arb_t x);

void _arb_get_mpfr(mpfr_t f, const fmpz_t mid, const fmpz_t exp, mpfr_rnd_t rnd);
void arb_get_mpfr(mpfr_t f, const arb_t x, mpfr_rnd_t rnd);

void arb_set_mpfr(arb_t y, const mpfr_t x, ulong error);

static __inline__ void
_arb_normalise(arb_t x)
{
    fmpz c = *arb_midref(x);

    if (COEFF_IS_MPZ(c))
    {
        long limbs, shift;
        __mpz_struct * z = COEFF_TO_PTR(c);

        limbs = mpz_size(z);
        shift = FLINT_BITS * (limbs-1) - arb_prec(x);

        if (shift > 0)
        {
            mpz_tdiv_q_2exp(z, z, shift);
            _fmpz_demote_val(arb_midref(x));
            fmpz_add_ui(arb_expref(x), arb_expref(x), shift);

            if (fmpz_is_zero(arb_radref(x)))
            {
                *arb_radref(x) = 1UL;
            }
            else
            {
                fmpz_cdiv_q_2exp(arb_radref(x), arb_radref(x), shift);
                fmpz_add_ui(arb_radref(x), arb_radref(x), 1UL);
            }
        }
    }
}

static __inline__ void
arb_zero(arb_t x)
{
    fmpz_zero(arb_midref(x));
    fmpz_zero(arb_expref(x));
    fmpz_zero(arb_radref(x));
}

void arb_set_si(arb_t x, long c);
void arb_set_ui(arb_t x, unsigned long c);
void arb_set_fmpz(arb_t x, const fmpz_t c);
void _arb_set_fmpq(arb_t x, const fmpz_t p, const fmpz_t q);
void arb_set_fmpq(arb_t x, const fmpq_t c);

void arb_add(arb_t z, const arb_t x, const arb_t y);
void arb_sub(arb_t z, const arb_t x, const arb_t y);
void arb_mul(arb_t z, const arb_t x, const arb_t y);
void arb_addmul(arb_t z, const arb_t x, const arb_t y);
void arb_submul(arb_t z, const arb_t x, const arb_t y);
void arb_div(arb_t c, const arb_t a, const arb_t b);

void arb_mul_si(arb_t y, const arb_t x, long c);

void arb_sqrt_fmpz(arb_t x, const fmpz_t n);
void arb_sqrt_ui(arb_t x, ulong n);

void arb_log_ui(arb_t x, ulong n);
void arb_log(arb_t y, const arb_t x);

void arb_exp(arb_t y, const arb_t x);

void arb_add_error_2exp(arb_t x, long c);
void _arb_rad_add_ufloat(arb_t y, const ufloat_t err);

void arb_const_pi_chudnovsky(arb_t x);
void arb_const_euler_brent_mcmillan(arb_t x);
void arb_const_zeta3_bsplit(arb_t x);

void arb_zeta_ui_bsplit(arb_t x, ulong s);
void arb_zeta_ui_mpfr(arb_t x, ulong n);
void arb_zeta_ui(arb_t x, ulong n);

/* fmpz extras */

static __inline__ void
_fmpz_set_si_small(fmpz_t x, long c)
{
    _fmpz_demote(x);
    *x = c;
}

static __inline__ void
_fmpz_mul_abs(fmpz_t f, const fmpz_t g, const fmpz_t h)
{
    fmpz_mul(f, g, h);
    fmpz_abs(f, f);
}

void _fmpz_addmul_abs(fmpz_t, const fmpz_t, const fmpz_t);

void _fmpz_addmul_abs_ui(fmpz_t, const fmpz_t, ulong x);

#endif
