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

#ifndef HYPGEOM_H
#define HYPGEOM_H

#include "fmprb.h"
#include "fmpz_poly.h"

typedef struct
{
    fmpz_poly_t A;
    fmpz_poly_t B;
    fmpz_poly_t P;
    fmpz_poly_t Q;

    /* precomputation data */
    int have_precomputed;
    long r;
    long boundC;
    long boundD;
    long boundK;
    fmpr_t MK;
}
hypgeom_struct;

typedef hypgeom_struct hypgeom_t[1];

void hypgeom_init(hypgeom_t hyp);

void hypgeom_clear(hypgeom_t hyp);

void hypgeom_precompute(hypgeom_t hyp);

long hypgeom_estimate_terms(const fmpr_t z, int r, long prec);

long hypgeom_bound(fmpr_t error, int r,
    long C, long D, long K, const fmpr_t TK, const fmpr_t z, long prec);

void fmprb_hypgeom_sum(fmprb_t P, fmprb_t Q, const hypgeom_t hyp, const long n, long prec);

void fmprb_hypgeom_infsum(fmprb_t P, fmprb_t Q, hypgeom_t hyp, long target_prec, long prec);

#endif

