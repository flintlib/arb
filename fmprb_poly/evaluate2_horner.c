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
_fmprb_poly_evaluate2_horner(fmprb_t y, fmprb_t z, fmprb_srcptr poly,
        long len, const fmprb_t x, long prec)
{
    if (len == 0)
    {
        fmprb_zero(y);
        fmprb_zero(z);
    }
    else if (len == 1)
    {
        fmprb_set_round(y, poly + 0, prec);
        fmprb_zero(z);
    }
    else if (fmprb_is_zero(x))
    {
        fmprb_set_round(y, poly + 0, prec);
        fmprb_set_round(z, poly + 1, prec);
    }
    else if (len == 2)
    {
        fmprb_mul(y, x, poly + 1, prec);
        fmprb_add(y, y, poly + 0, prec);
        fmprb_set_round(z, poly + 1, prec);
    }
    else
    {
        fmprb_t t, u, v;
        long i;

        fmprb_init(t);
        fmprb_init(u);
        fmprb_init(v);

        fmprb_set_round(u, poly + len - 1, prec);
        fmprb_zero(v);

        for (i = len - 2; i >= 0; i--)
        {
            fmprb_mul(t, v, x, prec);
            fmprb_add(v, u, t, prec);
            fmprb_mul(t, u, x, prec);
            fmprb_add(u, t, poly + i, prec);
        }

        fmprb_swap(y, u);
        fmprb_swap(z, v);

        fmprb_clear(t);
        fmprb_clear(u);
        fmprb_clear(v);
    }
}

void
fmprb_poly_evaluate2_horner(fmprb_t r, fmprb_t s, const fmprb_poly_t f, const fmprb_t a, long prec)
{
    _fmprb_poly_evaluate2_horner(r, s, f->coeffs, f->length, a, prec);
}

