/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#ifndef ARB_POLY_H
#define ARB_POLY_H

#include "arb.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"

static __inline__ void
fmpz_max(fmpz_t x, const fmpz_t a, const fmpz_t b)
{
    if (fmpz_cmp(a, b) >= 0)
        fmpz_set(x, a);
    else
        fmpz_set(x, b);
}

typedef struct
{
    fmpz * coeffs;
    long length;
    long alloc;
    fmpz exp;
    fmpz rad;
    long prec;
}
arb_poly_struct;

typedef arb_poly_struct arb_poly_t[1];

#define arb_poly_coeffs(x) ((x)->coeffs)
#define arb_poly_midref(x) (&(x)->mid)
#define arb_poly_expref(x) (&(x)->exp)
#define arb_poly_radref(x) (&(x)->rad)
#define arb_poly_prec(x) ((x)->prec)

void arb_poly_init(arb_poly_t poly, long prec);

void arb_poly_clear(arb_poly_t poly);

static __inline__ void
arb_poly_zero(arb_poly_t poly)
{
    poly->length = 0;         /* XXX: free memory? */
    fmpz_zero(arb_poly_radref(poly));
    fmpz_zero(arb_poly_expref(poly));
}

void arb_poly_set(arb_poly_t z, const arb_poly_t x);

void _arb_poly_fit_length(arb_poly_t poly, long length);

void _arb_poly_normalise(fmpz * h, fmpz_t hexp, fmpz_t hrad, long len, long bits);

static __inline__ void
arb_poly_set_si(arb_poly_t poly, long c)
{
    if (c == 0)
    {
        arb_poly_zero(poly);
    }
    else
    {
        _arb_poly_fit_length(poly, 1);
        fmpz_set_si(arb_poly_coeffs(poly), c);
        fmpz_zero(arb_poly_expref(poly));
        fmpz_zero(arb_poly_radref(poly));
        poly->length = 1;
    }
}

void arb_poly_debug(const arb_poly_t x);

void arb_poly_randtest(arb_poly_t x, flint_rand_t state, long len, long exp_bits);
void arb_poly_get_rand_fmpq_poly(fmpq_poly_t q, flint_rand_t state, const arb_poly_t x);

void arb_poly_set_fmpq_poly(arb_poly_t poly, const fmpq_poly_t src);

int arb_poly_contains_fmpq_poly(const arb_poly_t x, const fmpq_poly_t y);

void arb_poly_add(arb_poly_t z, const arb_poly_t x, const arb_poly_t y);
void arb_poly_neg(arb_poly_t z, const arb_poly_t x);
void arb_poly_mul(arb_poly_t z, const arb_poly_t x, const arb_poly_t y);
void arb_poly_mullow(arb_poly_t z, const arb_poly_t x, const arb_poly_t y, long n);

void arb_poly_inv_series(arb_poly_t z, const arb_poly_t x, long n);

void arb_poly_derivative(arb_poly_t z, const arb_poly_t x);
void arb_poly_integral(arb_poly_t z, const arb_poly_t x);

void arb_poly_log_series(arb_poly_t z, const arb_poly_t x, long n);
void arb_poly_exp_series(arb_poly_t z, const arb_poly_t x, long n);

#endif
