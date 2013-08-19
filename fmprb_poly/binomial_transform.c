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

    Copyright (C) 2013 Fredrik Johansson

******************************************************************************/

#include "fmprb_poly.h"

void
_fmprb_poly_binomial_transform(fmprb_ptr b, fmprb_srcptr a, long alen, long len, long prec)
{
    if (alen < 10 || len < 10)
        _fmprb_poly_binomial_transform_basecase(b, a, alen, len, prec);
    else
        _fmprb_poly_binomial_transform_convolution(b, a, alen, len, prec);
}

void
fmprb_poly_binomial_transform(fmprb_poly_t b, const fmprb_poly_t a, long len, long prec)
{
    if (len == 0 || a->length == 0)
    {
        fmprb_poly_zero(b);
        return;
    }

    if (b == a)
    {
        fmprb_poly_t c;
        fmprb_poly_init2(c, len);
        _fmprb_poly_binomial_transform(c->coeffs, a->coeffs, a->length, len, prec);
        fmprb_poly_swap(b, c);
        fmprb_poly_clear(c);
    }
    else
    {
        fmprb_poly_fit_length(b, len);
        _fmprb_poly_binomial_transform(b->coeffs, a->coeffs, a->length, len, prec);
    }

    _fmprb_poly_set_length(b, len);
    _fmprb_poly_normalise(b);
}

