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

#include "acb_poly.h"

void
acb_poly_set_arb_poly(acb_poly_t poly, const arb_poly_t re)
{
    long i, len;

    len = arb_poly_length(re);

    acb_poly_fit_length(poly, len);

    for (i = 0; i < len; i++)
    {
        arb_set(acb_realref(poly->coeffs + i), re->coeffs + i);
        arb_zero(acb_imagref(poly->coeffs + i));
    }

    _acb_poly_set_length(poly, len);
}

void
acb_poly_set2_arb_poly(acb_poly_t poly, const arb_poly_t re, const arb_poly_t im)
{
    long i, rlen, ilen, len;

    rlen = arb_poly_length(re);
    ilen = arb_poly_length(im);
    len = FLINT_MAX(rlen, ilen);

    acb_poly_fit_length(poly, len);

    for (i = 0; i < rlen; i++)
        arb_set(acb_realref(poly->coeffs + i), re->coeffs + i);
    for (i = rlen; i < len; i++)
        arb_zero(acb_realref(poly->coeffs + i));

    for (i = 0; i < ilen; i++)
        arb_set(acb_imagref(poly->coeffs + i), im->coeffs + i);
    for (i = ilen; i < len; i++)
        arb_zero(acb_imagref(poly->coeffs + i));

    _acb_poly_set_length(poly, len);
}
