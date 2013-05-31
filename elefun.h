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

    Copyright (C) 2013 Fredrik Johansson

******************************************************************************/

#ifndef ELEFUN_H
#define ELEFUN_H

#include "fmprb.h"
#include "fmprb_poly.h"
#include "fmpz_extras.h"

static __inline__ void
fmpz_print_fixed(const fmpz_t x, long prec)
{
    fmpr_t t;
    fmpr_init(t);
    fmpr_set_fmpz(t, x);
    fmpr_mul_2exp_si(t, t, -prec);
    fmpr_printd(t, prec / 3.33);
    fmpr_clear(t);
}

#define EXP_CACHE_PREC 1024
#define EXP_CACHE_BITS 8
#define EXP_CACHE_NUM (1L<<EXP_CACHE_BITS)
#define EXP_CACHE_LEVELS 2
#define EXP_CACHE_REDUCTION (EXP_CACHE_LEVELS * EXP_CACHE_BITS)

long elefun_exp_taylor_bound(long mag, long prec);

void elefun_exp_fixed_taylor_horner_precomp(fmpz_t y, fmpz_t yerr, const fmpz_t x, long n, long prec);

void elefun_exp_fixed_precomp(fmpz_t y, fmpz_t yerr, fmpz_t exponent,
    const fmpz_t x, const fmpz_t xerr, long prec);

int elefun_exp_precomp(fmprb_t z, const fmprb_t x, long prec, int minus_one);

void elefun_exp_via_mpfr(fmprb_t z, const fmprb_t x, long prec);


void _elefun_cos_minpoly_roots(fmprb_struct * alpha, long d, ulong n, long prec);

void _elefun_cos_minpoly(fmpz * coeffs, long d, ulong n);

void elefun_cos_minpoly(fmpz_poly_t poly, ulong n);

#endif

