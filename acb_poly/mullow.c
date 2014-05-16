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

#define CUTOFF 4

void
_acb_poly_mullow(acb_ptr res,
    acb_srcptr poly1, long len1,
    acb_srcptr poly2, long len2, long n, long prec)
{
    if (n < CUTOFF || len1 < CUTOFF || len2 < CUTOFF)
        _acb_poly_mullow_classical(res, poly1, len1, poly2, len2, n, prec);
    else
        _acb_poly_mullow_transpose(res, poly1, len1, poly2, len2, n, prec);
}

void
acb_poly_mullow(acb_poly_t res, const acb_poly_t poly1,
                                            const acb_poly_t poly2,
                                                long n, long prec)
{
    long len1, len2;

    len1 = poly1->length;
    len2 = poly2->length;

    if (len1 == 0 || len2 == 0 || n == 0)
    {
        acb_poly_zero(res);
        return;
    }

    n = FLINT_MIN((len1 + len2 - 1), n);
    len1 = FLINT_MIN(len1, n);
    len2 = FLINT_MIN(len2, n);

    if (res == poly1 || res == poly2)
    {
        acb_poly_t t;
        acb_poly_init2(t, n);
        _acb_poly_mullow(t->coeffs, poly1->coeffs, len1,
                                poly2->coeffs, len2, n, prec);
        acb_poly_swap(res, t);
        acb_poly_clear(t);
    }
    else
    {
        acb_poly_fit_length(res, n);
        _acb_poly_mullow(res->coeffs, poly1->coeffs, len1,
                                poly2->coeffs, len2, n, prec);
    }

    _acb_poly_set_length(res, n);
    _acb_poly_normalise(res);
}

