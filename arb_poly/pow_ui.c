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

#include "arb_poly.h"

void
_arb_poly_pow_ui(arb_ptr res, arb_srcptr f, long flen, ulong exp, long prec)
{
    _arb_poly_pow_ui_trunc_binexp(res, f, flen, exp, exp * (flen - 1) + 1, prec);
}

void
arb_poly_pow_ui(arb_poly_t res,
    const arb_poly_t poly, ulong exp, long prec)
{
    long flen, rlen;

    flen = poly->length;

    if (exp == 0)
    {
        arb_poly_one(res);
    }
    else if (flen == 0)
    {
        arb_poly_zero(res);
    }
    else
    {
        rlen = exp * (flen - 1) + 1;

        if (res != poly)
        {
            arb_poly_fit_length(res, rlen);
            _arb_poly_pow_ui(res->coeffs,
                poly->coeffs, flen, exp, prec);
            _arb_poly_set_length(res, rlen);
            _arb_poly_normalise(res);
        }
        else
        {
            arb_poly_t t;
            arb_poly_init2(t, rlen);
            _arb_poly_pow_ui(t->coeffs,
                poly->coeffs, flen, exp, prec);
            _arb_poly_set_length(t, rlen);
            _arb_poly_normalise(t);
            arb_poly_swap(res, t);
            arb_poly_clear(t);
        }
    }
}

