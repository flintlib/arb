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

#ifdef __cplusplus
}
#endif

#endif

