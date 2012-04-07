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

/* z = x + y/2^shift, assumes lengths > 0, shift >= 0 */
void
_arb_poly_add_shift(fmpz * z, fmpz_t zrad,
    const fmpz * x, long xlen, const fmpz_t xrad, 
    const fmpz * y, long ylen, const fmpz_t yrad, long shift)
{
    fmpz_t t;
    long i;

    fmpz_init(t);

    if (xlen >= ylen)
    {
        for (i = 0; i < ylen; i++)
        {
            fmpz_tdiv_q_2exp(t, y + i, shift);
            fmpz_add(z + i, x + i, t);
        }

        for (i = ylen; i < xlen; i++)
            fmpz_set(z + i, x + i);
    }
    else
    {
        for (i = 0; i < xlen; i++)
        {
            fmpz_tdiv_q_2exp(t, y + i, shift);
            fmpz_add(z + i, x + i, t);
        }

        for (i = xlen; i < ylen; i++)
            fmpz_tdiv_q_2exp(z + i, y + i, shift);
    }

    fmpz_cdiv_q_2exp(t, yrad, shift);
    fmpz_add(zrad, xrad, t);
    fmpz_add_ui(zrad, zrad, 1UL);

    fmpz_clear(t);
}


void
arb_poly_add(arb_poly_t z, const arb_poly_t x, const arb_poly_t y)
{
    long xlen, ylen, zlen;
    fmpz *zp, *xp, *yp;

    xlen = x->length;
    ylen = y->length;

    if (xlen == 0)
    {
        arb_poly_set(z, y);
        return;
    }

    if (ylen == 0)
    {
        arb_poly_set(z, x);
        return;
    }

    zlen = FLINT_MAX(xlen, ylen);
    _arb_poly_fit_length(z, zlen);

    zp = arb_poly_coeffs(z);
    xp = arb_poly_coeffs(x);
    yp = arb_poly_coeffs(y);

    if (fmpz_equal(arb_poly_expref(x), arb_poly_expref(y)))
    {
        _fmpz_poly_add(zp, xp, xlen, yp, ylen);
        fmpz_set(arb_poly_expref(z), arb_poly_expref(x));
        fmpz_add(arb_poly_radref(z), arb_poly_radref(x), arb_poly_radref(y));
    }
    else
    {
        /* TODO: handle huge exponents */
        /* TODO: be smart when one poly is exact */
        /* TODO: aliasing? */
        long shift;

        shift = fmpz_get_si(arb_poly_expref(x)) -
                    fmpz_get_si(arb_poly_expref(y));

        if (shift > 0)
        {
            fmpz_set(arb_poly_expref(z), arb_poly_expref(x));
            _arb_poly_add_shift(zp, arb_poly_radref(z),
                xp, xlen, arb_poly_radref(x),
                yp, ylen, arb_poly_radref(y), shift);
        }
        else
        {
            fmpz_set(arb_poly_expref(z), arb_poly_expref(y));
            _arb_poly_add_shift(zp, arb_poly_radref(z),
                yp, ylen, arb_poly_radref(y),
                xp, xlen, arb_poly_radref(x), -shift);
        }
    }

    z->length = zlen;
}
