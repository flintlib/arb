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

#ifndef ACB_HYPGEOM_H
#define ACB_HYPGEOM_H

#include "acb.h"
#include "acb_poly.h"

#ifdef __cplusplus
extern "C" {
#endif

void acb_hypgeom_pfq_bound_factor(mag_t C,
    acb_srcptr a, long p, acb_srcptr b, long q, const acb_t z, ulong n);

long acb_hypgeom_pfq_choose_n(acb_srcptr a, long p,
                         acb_srcptr b, long q, const acb_t z, long prec);

void acb_hypgeom_pfq_sum_forward(acb_t s, acb_t t, acb_srcptr a, long p, acb_srcptr b, long q,
    const acb_t z, long n, long prec);

void acb_hypgeom_pfq_sum_rs(acb_t s, acb_t t, acb_srcptr a, long p, acb_srcptr b, long q,
    const acb_t z, long n, long prec);

void acb_hypgeom_pfq_sum(acb_t s, acb_t t, acb_srcptr a, long p, acb_srcptr b, long q,
    const acb_t z, long n, long prec);

void acb_hypgeom_pfq_direct(acb_t res, acb_srcptr a, long p, acb_srcptr b, long q,
    const acb_t z, long n, long prec);

void acb_hypgeom_pfq_series_direct(acb_poly_t res,
    const acb_poly_struct * a, long p,
    const acb_poly_struct * b, long q,
    const acb_poly_t z,
    long n, long len, long prec, int regularized);

void acb_hypgeom_u_asymp(acb_t res, const acb_t a, const acb_t b,
    const acb_t z, long n, long prec);

void acb_hypgeom_erf_1f1a(acb_t res, const acb_t z, long prec);

void acb_hypgeom_erf_1f1b(acb_t res, const acb_t z, long prec);

void acb_hypgeom_erf_asymp(acb_t res, const acb_t z, long prec, long prec2);

void acb_hypgeom_erf(acb_t res, const acb_t z, long prec);

void acb_hypgeom_bessel_j_0f1(acb_t res, const acb_t nu, const acb_t z, long prec);

void acb_hypgeom_bessel_j_asymp(acb_t res, const acb_t nu, const acb_t z, long prec);

void acb_hypgeom_bessel_j(acb_t res, const acb_t nu, const acb_t z, long prec);

void acb_hypgeom_gamma_upper_asymp(acb_t res, const acb_t s, const acb_t z, int modified, long prec);

void acb_hypgeom_gamma_upper_1f1a(acb_t res, const acb_t s, const acb_t z, int modified, long prec);

void acb_hypgeom_gamma_upper_1f1b(acb_t res, const acb_t s, const acb_t z, int modified, long prec);

void acb_hypgeom_gamma_upper_singular(acb_t res, long s, const acb_t z, int modified, long prec);

void acb_hypgeom_gamma_upper(acb_t res, const acb_t s, const acb_t z, int modified, long prec);

void acb_hypgeom_expint(acb_t res, const acb_t s, const acb_t z, long prec);

void acb_hypgeom_erfc(acb_t res, const acb_t z, long prec);

void acb_hypgeom_erfi(acb_t res, const acb_t z, long prec);

#ifdef __cplusplus
}
#endif

#endif

