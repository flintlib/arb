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

#ifndef FMPCB_POLY_H
#define FMPCB_POLY_H

#include "fmpcb.h"
#include "fmprb_poly.h"

typedef struct
{
    fmpcb_struct * coeffs;
    long length;
    long alloc;
}
fmpcb_poly_struct;

typedef fmpcb_poly_struct fmpcb_poly_t[1];


/* Memory management */

void fmpcb_poly_init(fmpcb_poly_t poly);

void fmpcb_poly_init2(fmpcb_poly_t poly, long len);

void fmpcb_poly_clear(fmpcb_poly_t poly);

void fmpcb_poly_fit_length(fmpcb_poly_t poly, long len);

void _fmpcb_poly_set_length(fmpcb_poly_t poly, long len);

void _fmpcb_poly_normalise(fmpcb_poly_t poly);

static __inline__ void
fmpcb_poly_swap(fmpcb_poly_t poly1, fmpcb_poly_t poly2)
{
    fmpcb_poly_struct t = *poly1;
    *poly1 = *poly2;
    *poly2 = t;
}

static __inline__ long fmpcb_poly_length(const fmpcb_poly_t poly)
{
    return poly->length;
}

static __inline__ void fmpcb_poly_zero(fmpcb_poly_t poly)
{
    poly->length = 0;
}

static __inline__ void
fmpcb_poly_one(fmpcb_poly_t poly)
{
    fmpcb_poly_fit_length(poly, 1);
    fmpcb_one(poly->coeffs);
    _fmpcb_poly_set_length(poly, 1);
}

void fmpcb_poly_printd(const fmpcb_poly_t poly, long digits);

void _fmpcb_poly_evaluate(fmpcb_t res, const fmpcb_struct * f, long len, const fmpcb_t a, long prec);

void fmpcb_poly_evaluate(fmpcb_t res, const fmpcb_poly_t f, const fmpcb_t a, long prec);

void _fmpcb_poly_derivative(fmpcb_struct * res, const fmpcb_struct * poly, long len, long prec);

void fmpcb_poly_derivative(fmpcb_poly_t res, const fmpcb_poly_t poly, long prec);

void fmpcb_poly_set2_fmprb_poly(fmpcb_poly_t poly, const fmprb_poly_t re, const fmprb_poly_t im);

void fmpcb_poly_set2_fmpq_poly(fmpcb_poly_t poly, const fmpq_poly_t re, const fmpq_poly_t im, long prec);

void _fmpcb_poly_root_inclusion(fmpcb_t r, const fmpcb_t m,
    const fmpcb_struct * poly,
    const fmpcb_struct * polyder, long len, long prec);

#endif
