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
_arb_poly_evaluate2_rectangular(arb_t y, arb_t z, arb_srcptr poly,
    long len, const arb_t x, long prec)
{
    long i, j, m, r;
    arb_ptr xs;
    arb_t s, t, c;

    if (len < 3)
    {
        if (len == 0)
        {
            arb_zero(y);
            arb_zero(z);
        }
        else if (len == 1)
        {
            arb_set_round(y, poly + 0, prec);
            arb_zero(z);
        }
        else if (len == 2)
        {
            arb_mul(y, x, poly + 1, prec);
            arb_add(y, y, poly + 0, prec);
            arb_set_round(z, poly + 1, prec);
        }
        return;
    }

    m = n_sqrt(len) + 1;
    m *= 1;

    r = (len + m - 1) / m;

    xs = _arb_vec_init(m + 1);
    arb_init(s);
    arb_init(t);
    arb_init(c);

    _arb_vec_set_powers(xs, x, m + 1, prec);

    arb_set(y, poly + (r - 1) * m);
    for (j = 1; (r - 1) * m + j < len; j++)
        arb_addmul(y, xs + j, poly + (r - 1) * m + j, prec);

    for (i = r - 2; i >= 0; i--)
    {
        arb_set(s, poly + i * m);
        for (j = 1; j < m; j++)
            arb_addmul(s, xs + j, poly + i * m + j, prec);

        arb_mul(y, y, xs + m, prec);
        arb_add(y, y, s, prec);
    }

    len -= 1;
    r = (len + m - 1) / m;
    arb_mul_ui(z, poly + (r - 1) * m + 1, (r - 1) * m + 1, ARF_PREC_EXACT);
    for (j = 1; (r - 1) * m + j < len; j++)
    {
        arb_mul_ui(c, poly + (r - 1) * m + j + 1, (r - 1) * m + j + 1, ARF_PREC_EXACT);
        arb_addmul(z, xs + j, c, prec);
    }

    for (i = r - 2; i >= 0; i--)
    {
        arb_mul_ui(s, poly + i * m + 1, i * m + 1, ARF_PREC_EXACT);

        for (j = 1; j < m; j++)
        {
            arb_mul_ui(c, poly + i * m + j + 1, i * m + j + 1, ARF_PREC_EXACT);
            arb_addmul(s, xs + j, c, prec);
        }

        arb_mul(z, z, xs + m, prec);
        arb_add(z, z, s, prec);
    }

    _arb_vec_clear(xs, m + 1);
    arb_clear(s);
    arb_clear(t);
    arb_clear(c);
}

void
arb_poly_evaluate2_rectangular(arb_t r, arb_t s, const arb_poly_t f, const arb_t a, long prec)
{
    _arb_poly_evaluate2_rectangular(r, s, f->coeffs, f->length, a, prec);
}

