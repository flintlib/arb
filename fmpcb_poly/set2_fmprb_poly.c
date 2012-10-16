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

#include "fmpcb_poly.h"

void
fmpcb_poly_set_fmprb_poly(fmpcb_poly_t poly, const fmprb_poly_t re)
{
    long i, len;

    len = fmprb_poly_length(re);

    fmpcb_poly_fit_length(poly, len);

    for (i = 0; i < len; i++)
    {
        fmprb_set(fmpcb_realref(poly->coeffs + i), re->coeffs + i);
        fmprb_zero(fmpcb_imagref(poly->coeffs + i));
    }

    _fmpcb_poly_set_length(poly, len);
}

void
fmpcb_poly_set2_fmprb_poly(fmpcb_poly_t poly, const fmprb_poly_t re, const fmprb_poly_t im)
{
    long i, rlen, ilen, len;

    rlen = fmprb_poly_length(re);
    ilen = fmprb_poly_length(im);
    len = FLINT_MAX(rlen, ilen);

    fmpcb_poly_fit_length(poly, len);

    for (i = 0; i < rlen; i++)
        fmprb_set(fmpcb_realref(poly->coeffs + i), re->coeffs + i);
    for (i = rlen; i < len; i++)
        fmprb_zero(fmpcb_realref(poly->coeffs + i));

    for (i = 0; i < ilen; i++)
        fmprb_set(fmpcb_imagref(poly->coeffs + i), im->coeffs + i);
    for (i = ilen; i < len; i++)
        fmprb_zero(fmpcb_imagref(poly->coeffs + i));

    _fmpcb_poly_set_length(poly, len);
}
