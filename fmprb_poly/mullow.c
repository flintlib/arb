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

#include "fmprb_poly.h"

#define BLOCK_CUTOFF 5
#define SCALE_CUTOFF 50

#define BLOCK2_CUTOFF 7

void
_fmprb_poly_mullow(fmprb_ptr res,
    fmprb_srcptr poly1, long len1,
    fmprb_srcptr poly2, long len2, long n, long prec)
{
    if (n == 1)
    {
        fmprb_mul(res, poly1, poly2, prec);
    }
    else
    {
#if 0
        if (n < BLOCK_CUTOFF || len1 < BLOCK_CUTOFF || len2 < BLOCK_CUTOFF)
            _fmprb_poly_mullow_classical(res, poly1, len1, poly2, len2, n, prec);
        else if (n < SCALE_CUTOFF || len1 < SCALE_CUTOFF || len2 < SCALE_CUTOFF)
            _fmprb_poly_mullow_block(res, poly1, len1, poly2, len2, n, prec);
        else
            _fmprb_poly_mullow_block_scaled(res, poly1, len1, poly2, len2, n, prec);
#else
        if (n < BLOCK2_CUTOFF || len1 < BLOCK2_CUTOFF || len2 < BLOCK2_CUTOFF)
            _fmprb_poly_mullow_classical(res, poly1, len1, poly2, len2, n, prec);
        else
            _fmprb_poly_mullow_block2(res, poly1, len1, poly2, len2, n, prec);
#endif
    }
}

void
fmprb_poly_mullow(fmprb_poly_t res, const fmprb_poly_t poly1,
                                            const fmprb_poly_t poly2,
                                                long n, long prec)
{
    long len_out;

    if (poly1->length == 0 || poly2->length == 0 || n == 0)
    {
        fmprb_poly_zero(res);
        return;
    }

    len_out = poly1->length + poly2->length - 1;
    if (n > len_out)
        n = len_out;

    if (res == poly1 || res == poly2)
    {
        fmprb_poly_t t;
        fmprb_poly_init2(t, n);
        _fmprb_poly_mullow(t->coeffs, poly1->coeffs, poly1->length,
                                poly2->coeffs, poly2->length, n, prec);
        fmprb_poly_swap(res, t);
        fmprb_poly_clear(t);
    }
    else
    {
        fmprb_poly_fit_length(res, n);
        _fmprb_poly_mullow(res->coeffs, poly1->coeffs, poly1->length,
                                poly2->coeffs, poly2->length, n, prec);
    }

    _fmprb_poly_set_length(res, n);
    _fmprb_poly_normalise(res);
}
