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
_arb_poly_mul(fmpz * h, fmpz_t hexp, fmpz_t hrad,
            const fmpz * f, const fmpz_t fexp, const fmpz_t frad, long flen,
            const fmpz * g, const fmpz_t gexp, const fmpz_t grad, long glen,
            long maxbits)
{
    /* Multiply errors */
    if (fmpz_is_zero(frad) && fmpz_is_zero(grad))
    {
        fmpz_zero(hrad);
    }
    else
    {
        long terms, fbits, gbits;
        fmpz_t t, u;

        terms = FLINT_MIN(flen, glen);
        fbits = _fmpz_vec_max_bits(f, flen);
        gbits = _fmpz_vec_max_bits(g, glen);
        fbits = FLINT_ABS(fbits);
        gbits = FLINT_ABS(gbits);

        fmpz_init(t);
        fmpz_init(u);

        /* For a single coefficient product, we have
           (a + r)(b + s) = a*b + s*a + r*b + r*s */
        fmpz_mul_2exp(t, frad, gbits);  /* r*b */
        fmpz_mul_2exp(u, grad, fbits);  /* s*a */
        fmpz_add(t, t, u);
        fmpz_addmul(t, frad, grad);     /* r*s */

        /* add up for all terms */
        fmpz_mul_ui(hrad, t, terms);

        fmpz_clear(t);
        fmpz_clear(u);
    }

    /* Exact multiplication */
    if (flen >= glen)
        _fmpz_poly_mul(h, f, flen, g, glen);
    else
        _fmpz_poly_mul(h, g, glen, f, flen);

    fmpz_add(hexp, fexp, gexp);

    _arb_poly_normalise(h, hexp, hrad, flen + glen - 1, maxbits);
}

void
arb_poly_mul(arb_poly_t z, const arb_poly_t x, const arb_poly_t y)
{
    long xlen, ylen;

    xlen = x->length;
    ylen = y->length;

    if (xlen == 0 || ylen == 0)
    {
        arb_poly_zero(z);
        return;
    }

    if (z == x)
    {
        arb_poly_t t;
        arb_poly_init(t, arb_poly_prec(x));
        arb_poly_set(t, x);
        arb_poly_mul(z, t, y);
        arb_poly_clear(t);
        return;
    }

    if (z == y)
    {
        arb_poly_t t;
        arb_poly_init(t, arb_poly_prec(y));
        arb_poly_set(t, y);
        arb_poly_mul(z, x, t);
        arb_poly_clear(t);
        return;
    }

    _arb_poly_fit_length(z, xlen + ylen - 1);
    _arb_poly_mul(z->coeffs, &z->exp, &z->rad,
        x->coeffs, &x->exp, &x->rad, xlen,
        y->coeffs, &y->exp, &y->rad, ylen,
        z->prec);
    z->length = xlen + ylen - 1;
}
