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

#include "arb_poly.h"

#define BLOCK_CUTOFF 16

void
_arb_poly_mullow(arb_ptr res,
    arb_srcptr poly1, long len1,
    arb_srcptr poly2, long len2, long n, long prec)
{
    if (n == 1)
    {
        arb_mul(res, poly1, poly2, prec);
    }
    else
    {
        if (n < BLOCK_CUTOFF || len1 < BLOCK_CUTOFF || len2 < BLOCK_CUTOFF)
            _arb_poly_mullow_classical(res, poly1, len1, poly2, len2, n, prec);
        else
            _arb_poly_mullow_block(res, poly1, len1, poly2, len2, n, prec);
    }
}

void
arb_poly_mullow(arb_poly_t res, const arb_poly_t poly1,
                                            const arb_poly_t poly2,
                                                long n, long prec)
{
    long len_out;

    if (poly1->length == 0 || poly2->length == 0 || n == 0)
    {
        arb_poly_zero(res);
        return;
    }

    len_out = poly1->length + poly2->length - 1;
    if (n > len_out)
        n = len_out;

    if (res == poly1 || res == poly2)
    {
        arb_poly_t t;
        arb_poly_init2(t, n);
        _arb_poly_mullow(t->coeffs, poly1->coeffs, poly1->length,
                                poly2->coeffs, poly2->length, n, prec);
        arb_poly_swap(res, t);
        arb_poly_clear(t);
    }
    else
    {
        arb_poly_fit_length(res, n);
        _arb_poly_mullow(res->coeffs, poly1->coeffs, poly1->length,
                                poly2->coeffs, poly2->length, n, prec);
    }

    _arb_poly_set_length(res, n);
    _arb_poly_normalise(res);
}
