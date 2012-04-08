/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "arb_poly.h"

void
_arb_poly_mullow(fmpz * h, fmpz_t hexp, fmpz_t hrad,
            const fmpz * f, const fmpz_t fexp, const fmpz_t frad, long flen,
            const fmpz * g, const fmpz_t gexp, const fmpz_t grad, long glen,
            long maxbits, long trunc)
{
    /* Multiply errors */
    if (fmpz_is_zero(frad) && fmpz_is_zero(grad))
    {
        fmpz_zero(hrad);
    }
    else
    {
        long terms;
        fmpz_t t;

        terms = FLINT_MIN(flen, glen);
        fmpz_init(t);

        /* For a single coefficient product, we have
           (a + r)(b + s) = a*b + s*a + r*b + r*s */
        fmpz_mul(t, frad, grad);
        _fmpz_addmul_abs(t, frad, g + _fmpz_vec_height_index(g, glen));
        _fmpz_addmul_abs(t, grad, f + _fmpz_vec_height_index(f, flen));
        /* add up for all terms */
        fmpz_mul_ui(hrad, t, terms);

        fmpz_clear(t);
    }

    /* Exact multiplication */
    if (flen >= glen)
        _fmpz_poly_mullow(h, f, flen, g, glen, trunc);
    else
        _fmpz_poly_mullow(h, g, glen, f, flen, trunc);

    fmpz_add(hexp, fexp, gexp);

    _arb_poly_normalise(h, hexp, hrad, trunc, maxbits);
}

void
arb_poly_mullow(arb_poly_t z, const arb_poly_t x, const arb_poly_t y,
    long n)
{
    long xlen, ylen;

    xlen = x->length;
    ylen = y->length;

    if (xlen == 0 || ylen == 0 || n == 0)
    {
        arb_poly_zero(z);
        return;
    }

    if (z == x)
    {
        arb_poly_t t;
        arb_poly_init(t, arb_poly_prec(x));
        arb_poly_set(t, x);
        arb_poly_mullow(z, t, y, n);
        arb_poly_clear(t);
        return;
    }

    if (z == y)
    {
        arb_poly_t t;
        arb_poly_init(t, arb_poly_prec(y));
        arb_poly_set(t, y);
        arb_poly_mullow(z, x, t, n);
        arb_poly_clear(t);
        return;
    }

    xlen = FLINT_MIN(xlen, n);
    ylen = FLINT_MIN(ylen, n);
    n = FLINT_MIN(n, xlen + ylen - 1);

    _arb_poly_fit_length(z, n);
    _arb_poly_mullow(z->coeffs, &z->exp, &z->rad,
        x->coeffs, &x->exp, &x->rad, xlen,
        y->coeffs, &y->exp, &y->rad, ylen,
        z->prec, n);
    z->length = n;
}
