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
_fmprb_poly_evaluate2_fmpcb_horner(fmpcb_t y, fmpcb_t z,
    fmprb_srcptr poly, long len, const fmpcb_t x, long prec)
{
    if (len == 0)
    {
        fmpcb_zero(y);
        fmpcb_zero(z);
    }
    else if (len == 1)
    {
        fmpcb_set_round_fmprb(y, poly + 0, prec);
        fmpcb_zero(z);
    }
    else if (fmpcb_is_zero(x))
    {
        fmpcb_set_round_fmprb(y, poly + 0, prec);
        fmpcb_set_round_fmprb(z, poly + 1, prec);
    }
    else if (len == 2)
    {
        fmpcb_mul_fmprb(y, x, poly + 1, prec);
        fmpcb_add_fmprb(y, y, poly + 0, prec);
        fmpcb_set_round_fmprb(z, poly + 1, prec);
    }
    else
    {
        fmpcb_t t, u, v;
        long i;

        fmpcb_init(t);
        fmpcb_init(u);
        fmpcb_init(v);

        fmpcb_set_round_fmprb(u, poly + len - 1, prec);
        fmpcb_zero(v);

        for (i = len - 2; i >= 0; i--)
        {
            fmpcb_mul(t, v, x, prec);
            fmpcb_add(v, u, t, prec);
            fmpcb_mul(t, u, x, prec);
            fmpcb_add_fmprb(u, t, poly + i, prec);
        }

        fmpcb_swap(y, u);
        fmpcb_swap(z, v);

        fmpcb_clear(t);
        fmpcb_clear(u);
        fmpcb_clear(v);
    }
}

void
fmprb_poly_evaluate2_fmpcb_horner(fmpcb_t r, fmpcb_t s, const fmprb_poly_t f, const fmpcb_t a, long prec)
{
    _fmprb_poly_evaluate2_fmpcb_horner(r, s, f->coeffs, f->length, a, prec);
}

