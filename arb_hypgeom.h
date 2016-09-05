/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef ARB_HYPGEOM_H
#define ARB_HYPGEOM_H

#include "arb.h"
#include "arb_poly.h"

#ifdef __cplusplus
extern "C" {
#endif

void arb_hypgeom_pfq(arb_t res, arb_srcptr a, slong p, arb_srcptr b, slong q,
    const arb_t z, int regularized, slong prec);

void arb_hypgeom_0f1(arb_t res, const arb_t a, const arb_t z, int regularized, slong prec);
void arb_hypgeom_m(arb_t res, const arb_t a, const arb_t b, const arb_t z, int regularized, slong prec);
void arb_hypgeom_1f1(arb_t res, const arb_t a, const arb_t b, const arb_t z, int regularized, slong prec);
void arb_hypgeom_u(arb_t res, const arb_t a, const arb_t b, const arb_t z, slong prec);
void arb_hypgeom_2f1(arb_t res, const arb_t a, const arb_t b, const arb_t c, const arb_t z, int regularized, slong prec);

void arb_hypgeom_erf(arb_t res, const arb_t z, slong prec);
void _arb_hypgeom_erf_series(arb_ptr g, arb_srcptr h, slong hlen, slong len, slong prec);
void arb_hypgeom_erf_series(arb_poly_t g, const arb_poly_t h, slong len, slong prec);

void arb_hypgeom_erfc(arb_t res, const arb_t z, slong prec);
void _arb_hypgeom_erfc_series(arb_ptr g, arb_srcptr h, slong hlen, slong len, slong prec);
void arb_hypgeom_erfc_series(arb_poly_t g, const arb_poly_t h, slong len, slong prec);

void arb_hypgeom_erfi(arb_t res, const arb_t z, slong prec);
void _arb_hypgeom_erfi_series(arb_ptr g, arb_srcptr h, slong hlen, slong len, slong prec);
void arb_hypgeom_erfi_series(arb_poly_t g, const arb_poly_t h, slong len, slong prec);

void arb_hypgeom_fresnel(arb_t res1, arb_t res2, const arb_t z, int normalized, slong prec);
void _arb_hypgeom_fresnel_series(arb_ptr s, arb_ptr c, arb_srcptr h, slong hlen, int normalized, slong len, slong prec);
void arb_hypgeom_fresnel_series(arb_poly_t s, arb_poly_t c, const arb_poly_t h, int normalized, slong len, slong prec);

void arb_hypgeom_ei(arb_t res, const arb_t z, slong prec);
void _arb_hypgeom_ei_series(arb_ptr g, arb_srcptr h, slong hlen, slong len, slong prec);
void arb_hypgeom_ei_series(arb_poly_t g, const arb_poly_t h, slong len, slong prec);

void arb_hypgeom_si(arb_t res, const arb_t z, slong prec);
void _arb_hypgeom_si_series(arb_ptr g, arb_srcptr h, slong hlen, slong len, slong prec);
void arb_hypgeom_si_series(arb_poly_t g, const arb_poly_t h, slong len, slong prec);

void arb_hypgeom_ci(arb_t res, const arb_t z, slong prec);
void _arb_hypgeom_ci_series(arb_ptr g, arb_srcptr h, slong hlen, slong len, slong prec);
void arb_hypgeom_ci_series(arb_poly_t g, const arb_poly_t h, slong len, slong prec);

void arb_hypgeom_shi(arb_t res, const arb_t z, slong prec);
void _arb_hypgeom_shi_series(arb_ptr g, arb_srcptr h, slong hlen, slong len, slong prec);
void arb_hypgeom_shi_series(arb_poly_t g, const arb_poly_t h, slong len, slong prec);

void arb_hypgeom_chi(arb_t res, const arb_t z, slong prec);
void _arb_hypgeom_chi_series(arb_ptr g, arb_srcptr h, slong hlen, slong len, slong prec);
void arb_hypgeom_chi_series(arb_poly_t g, const arb_poly_t h, slong len, slong prec);

void arb_hypgeom_li(arb_t res, const arb_t z, int offset, slong prec);
void _arb_hypgeom_li_series(arb_ptr g, arb_srcptr h, slong hlen, int offset, slong len, slong prec);
void arb_hypgeom_li_series(arb_poly_t g, const arb_poly_t h, int offset, slong len, slong prec);

void arb_hypgeom_airy(arb_t ai, arb_t aip, arb_t bi, arb_t bip, const arb_t z, slong prec);
void arb_hypgeom_airy_jet(arb_ptr ai, arb_ptr bi, const arb_t z, slong len, slong prec);
void arb_hypgeom_airy_series(arb_poly_t ai, arb_poly_t ai_prime,
    arb_poly_t bi, arb_poly_t bi_prime, const arb_poly_t z, slong len, slong prec);
void _arb_hypgeom_airy_series(arb_ptr ai, arb_ptr ai_prime,
    arb_ptr bi, arb_ptr bi_prime, arb_srcptr z, slong zlen, slong len, slong prec);

#ifdef __cplusplus
}
#endif

#endif

