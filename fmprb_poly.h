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

#ifndef FMPRB_POLY_H
#define FMPRB_POLY_H

#include "fmprb.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"

typedef struct
{
    fmprb_struct * coeffs;
    long length;
    long alloc;
}
fmprb_poly_struct;

typedef fmprb_poly_struct fmprb_poly_t[1];


/* Memory management */

void fmprb_poly_init(fmprb_poly_t poly);

void fmprb_poly_clear(fmprb_poly_t poly);

void fmprb_poly_fit_length(fmprb_poly_t poly, long len);

void _fmprb_poly_set_length(fmprb_poly_t poly, long len);

void _fmprb_poly_normalise(fmprb_poly_t poly);

/* Basic manipulation */

static __inline__ void fmprb_poly_zero(fmprb_poly_t poly)
{
    poly->length = 0;
}

/* Conversions */

void fmprb_poly_set_fmpz_poly(fmprb_poly_t poly, const fmpz_poly_t src, long prec);

void fmprb_poly_set_fmpq_poly(fmprb_poly_t poly, const fmpq_poly_t src, long prec);

/* Comparisons */

int fmprb_poly_contains_fmpq_poly(const fmprb_poly_t poly1, const fmpq_poly_t poly2);

/* IO */

void fmprb_poly_printd(const fmprb_poly_t poly, long digits);

/* Arithmetic */

void fmprb_poly_add(fmprb_poly_t res, const fmprb_poly_t poly1,
              const fmprb_poly_t poly2, long prec);

void _fmprb_poly_mullow(fmprb_struct * C,
    const fmprb_struct * A, long lenA,
    const fmprb_struct * B, long lenB, long n, long prec);

void fmprb_poly_mullow(fmprb_poly_t res, const fmprb_poly_t poly1,
              const fmprb_poly_t poly2, long len, long prec);

void _fmprb_poly_mul(fmprb_struct * C,
    const fmprb_struct * A, long lenA,
    const fmprb_struct * B, long lenB, long prec);

void fmprb_poly_mul(fmprb_poly_t res, const fmprb_poly_t poly1,
              const fmprb_poly_t poly2, long prec);

void _fmprb_poly_inv_series(fmprb_struct * Qinv, const fmprb_struct * Q, long len, long prec);

void fmprb_poly_inv_series(fmprb_poly_t Qinv, const fmprb_poly_t Q, long n, long prec);

/* Special functions */

void fmprb_poly_log_gamma_series(fmprb_poly_t z, long n, long prec);

#endif

