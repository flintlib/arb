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

long _fmprb_bits(const fmprb_t x);

long _fmprb_vec_bits(const fmprb_struct * x, long len);

void
_fmprb_poly_evaluate2(fmprb_t y, fmprb_t z, const fmprb_struct * f, long len, const fmprb_t x, long prec)
{
    if ((prec >= 1024) && (len >= 5 + 20000 / prec))
    {
        long fbits, xbits;

        xbits = _fmprb_bits(x);
        fbits = _fmprb_vec_bits(f, len);

        if (fbits <= prec / 2)
        {
            _fmprb_poly_evaluate2_rectangular(y, z, f, len, x, prec);
            return;
        }
    }

    _fmprb_poly_evaluate2_horner(y, z, f, len, x, prec);
}

void
fmprb_poly_evaluate2(fmprb_t r, fmprb_t s, const fmprb_poly_t f, const fmprb_t a, long prec)
{
    _fmprb_poly_evaluate2(r, s, f->coeffs, f->length, a, prec);
}

