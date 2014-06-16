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

#include "arb.h"

void arb_pow_fmpz_binexp(arb_t y, const arb_t b, const fmpz_t e, long prec);
void arb_pow_fmpz(arb_t y, const arb_t b, const fmpz_t e, long prec);
void arb_pow_ui(arb_t y, const arb_t b, ulong e, long prec);
void arb_ui_pow_ui(arb_t y, ulong b, ulong e, long prec);
void arb_si_pow_ui(arb_t y, long b, ulong e, long prec);
void arb_pow_fmpq(arb_t y, const arb_t x, const fmpq_t a, long prec);

void arb_div_2expm1_ui(arb_t z, const arb_t x, ulong n, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    arb_get_fmprb(t, x);
    fmprb_div_2expm1_ui(t, t, n, prec);
    arb_set_fmprb(z, t);
    fmprb_clear(t);
}

void arb_pow(arb_t z, const arb_t x, const arb_t y, long prec)
{
    fmprb_t t, u;
    fmprb_init(t);
    fmprb_init(u);
    arb_get_fmprb(t, x);
    arb_get_fmprb(u, y);
    fmprb_pow(t, t, u, prec);
    arb_set_fmprb(z, t);
    fmprb_clear(t);
    fmprb_clear(u);
}

void arb_root(arb_t z, const arb_t x, ulong k, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    arb_get_fmprb(t, x);
    fmprb_root(t, t, k, prec);
    arb_set_fmprb(z, t);
    fmprb_clear(t);
}

void arb_log(arb_t z, const arb_t x, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    arb_get_fmprb(t, x);
    fmprb_log(t, t, prec);
    arb_set_fmprb(z, t);
    fmprb_clear(t);
}

void arb_log_ui(arb_t z, ulong x, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    fmprb_log_ui(t, x, prec);
    arb_set_fmprb(z, t);
    fmprb_clear(t);
}

void arb_log_fmpz(arb_t z, const fmpz_t x, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    fmprb_log_fmpz(t, x, prec);
    arb_set_fmprb(z, t);
    fmprb_clear(t);
}

void arb_exp(arb_t z, const arb_t x, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    arb_get_fmprb(t, x);
    fmprb_exp(t, t, prec);
    arb_set_fmprb(z, t);
    fmprb_clear(t);
}

void arb_expm1(arb_t z, const arb_t x, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    arb_get_fmprb(t, x);
    fmprb_expm1(t, t, prec);
    arb_set_fmprb(z, t);
    fmprb_clear(t);
}

void arb_sin_pi(arb_t z, const arb_t x, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    arb_get_fmprb(t, x);
    fmprb_sin_pi(t, t, prec);
    arb_set_fmprb(z, t);
    fmprb_clear(t);
}

void arb_cos_pi(arb_t z, const arb_t x, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    arb_get_fmprb(t, x);
    fmprb_cos_pi(t, t, prec);
    arb_set_fmprb(z, t);
    fmprb_clear(t);
}

void arb_sin_cos_pi(arb_t s, arb_t c, const arb_t x, long prec)
{
    fmprb_t t, u;
    fmprb_init(t);
    fmprb_init(u);
    arb_get_fmprb(t, x);
    fmprb_sin_cos_pi(t, u, t, prec);
    arb_set_fmprb(s, t);
    arb_set_fmprb(c, u);
    fmprb_clear(t);
    fmprb_clear(u);
}

void arb_sin_cos_pi_fmpq(arb_t s, arb_t c, const fmpq_t x, long prec)
{
    fmprb_t t, u;
    fmprb_init(t);
    fmprb_init(u);
    fmprb_sin_cos_pi_fmpq(t, u, x, prec);
    arb_set_fmprb(s, t);
    arb_set_fmprb(c, u);
    fmprb_clear(t);
    fmprb_clear(u);
}

void arb_sin_pi_fmpq(arb_t z, const fmpq_t x, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    fmprb_sin_pi_fmpq(t, x, prec);
    arb_set_fmprb(z, t);
    fmprb_clear(t);
}

void arb_cos_pi_fmpq(arb_t z, const fmpq_t x, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    fmprb_cos_pi_fmpq(t, x, prec);
    arb_set_fmprb(z, t);
    fmprb_clear(t);
}

void arb_atan(arb_t z, const arb_t x, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    arb_get_fmprb(t, x);
    fmprb_atan(t, t, prec);
    arb_set_fmprb(z, t);
    fmprb_clear(t);
}

void arb_atan2(arb_t z, const arb_t b, const arb_t a, long prec)
{
    fmprb_t t, u;
    fmprb_init(t);
    fmprb_init(u);
    arb_get_fmprb(t, b);
    arb_get_fmprb(u, a);
    fmprb_atan2(t, t, u, prec);
    arb_set_fmprb(z, t);
    fmprb_clear(t);
    fmprb_clear(u);
}

void arb_fac_ui(arb_t z, ulong n, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    fmprb_fac_ui(t, n, prec);
    arb_set_fmprb(z, t);
    fmprb_clear(t);
}

void arb_bin_ui(arb_t z, const arb_t n, ulong k, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    arb_get_fmprb(t, n);
    fmprb_bin_ui(t, t, k, prec);
    arb_set_fmprb(z, t);
    fmprb_clear(t);
}

void arb_bin_uiui(arb_t z, ulong n, ulong k, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    fmprb_bin_uiui(t, n, k, prec);
    arb_set_fmprb(z, t);
    fmprb_clear(t);
}

void arb_fib_fmpz(arb_t z, const fmpz_t n, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    fmprb_fib_fmpz(t, n, prec);
    arb_set_fmprb(z, t);
    fmprb_clear(t);
}

void arb_fib_ui(arb_t z, ulong n, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    fmprb_fib_ui(t, n, prec);
    arb_set_fmprb(z, t);
    fmprb_clear(t);
}

