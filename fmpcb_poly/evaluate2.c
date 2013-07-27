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

#include "fmpcb_poly.h"

void
_fmpcb_poly_evaluate2(fmpcb_t y, fmpcb_t z, fmpcb_srcptr f, long len, const fmpcb_t x, long prec)
{
    if ((prec >= 1024) && (len >= 5 + 20000 / prec))
    {
        long fbits;

        fbits = _fmpcb_vec_bits(f, len);

        if (fbits <= prec / 2)
        {
            _fmpcb_poly_evaluate2_rectangular(y, z, f, len, x, prec);
            return;
        }
    }

    _fmpcb_poly_evaluate2_horner(y, z, f, len, x, prec);
}

void
fmpcb_poly_evaluate2(fmpcb_t r, fmpcb_t s, const fmpcb_poly_t f, const fmpcb_t a, long prec)
{
    _fmpcb_poly_evaluate2(r, s, f->coeffs, f->length, a, prec);
}

