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

    Copyright (C) 2011 William Hart
    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "fmpcb_poly.h"

void
_fmpcb_poly_product_roots(fmpcb_ptr poly, fmpcb_srcptr xs, long n, long prec)
{
    if (n == 0)
    {
        fmpcb_one(poly);
    }
    else if (n == 1)
    {
        fmpcb_neg(poly, xs);
        fmpcb_one(poly + 1);
    }
    else if (n == 2)
    {
        fmpcb_mul(poly, xs + 0, xs + 1, prec);
        fmpcb_add(poly + 1, xs + 0, xs + 1, prec);
        fmpcb_neg(poly + 1, poly + 1);
        fmpcb_one(poly + 2);
    }
    else
    {
        const long m = (n + 1) / 2;
        fmpcb_ptr tmp;

        tmp = _fmpcb_vec_init(n + 2);

        _fmpcb_poly_product_roots(tmp, xs, m, prec);
        _fmpcb_poly_product_roots(tmp + m + 1, xs + m, n - m, prec);
        _fmpcb_poly_mul_monic(poly, tmp, m + 1, tmp + m + 1, n - m + 1, prec);

        _fmpcb_vec_clear(tmp, n + 2);
    }
}

void
fmpcb_poly_product_roots(fmpcb_poly_t poly, fmpcb_srcptr xs, long n, long prec)
{
    fmpcb_poly_fit_length(poly, n + 1);
    _fmpcb_poly_product_roots(poly->coeffs, xs, n, prec);
    _fmpcb_poly_set_length(poly, n + 1);
}
