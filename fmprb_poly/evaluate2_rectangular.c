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
_fmprb_poly_evaluate2_rectangular(fmprb_t y, fmprb_t z, fmprb_srcptr poly,
    long len, const fmprb_t x, long prec)
{
    long i, j, m, r;
    fmprb_ptr xs;
    fmprb_t s, t, c;

    if (len < 3)
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
        else if (len == 2)
        {
            fmprb_mul(y, x, poly + 1, prec);
            fmprb_add(y, y, poly + 0, prec);
            fmprb_set_round(z, poly + 1, prec);
        }
        return;
    }

    m = n_sqrt(len) + 1;
    m *= 1;

    r = (len + m - 1) / m;

    xs = _fmprb_vec_init(m + 1);
    fmprb_init(s);
    fmprb_init(t);
    fmprb_init(c);

    _fmprb_vec_set_powers(xs, x, m + 1, prec);

    fmprb_set(y, poly + (r - 1) * m);
    for (j = 1; (r - 1) * m + j < len; j++)
        fmprb_addmul(y, xs + j, poly + (r - 1) * m + j, prec);

    for (i = r - 2; i >= 0; i--)
    {
        fmprb_set(s, poly + i * m);
        for (j = 1; j < m; j++)
            fmprb_addmul(s, xs + j, poly + i * m + j, prec);

        fmprb_mul(y, y, xs + m, prec);
        fmprb_add(y, y, s, prec);
    }

    len -= 1;
    r = (len + m - 1) / m;
    fmprb_mul_ui(z, poly + (r - 1) * m + 1, (r - 1) * m + 1, FMPR_PREC_EXACT);
    for (j = 1; (r - 1) * m + j < len; j++)
    {
        fmprb_mul_ui(c, poly + (r - 1) * m + j + 1, (r - 1) * m + j + 1, FMPR_PREC_EXACT);
        fmprb_addmul(z, xs + j, c, prec);
    }

    for (i = r - 2; i >= 0; i--)
    {
        fmprb_mul_ui(s, poly + i * m + 1, i * m + 1, FMPR_PREC_EXACT);

        for (j = 1; j < m; j++)
        {
            fmprb_mul_ui(c, poly + i * m + j + 1, i * m + j + 1, FMPR_PREC_EXACT);
            fmprb_addmul(s, xs + j, c, prec);
        }

        fmprb_mul(z, z, xs + m, prec);
        fmprb_add(z, z, s, prec);
    }

    _fmprb_vec_clear(xs, m + 1);
    fmprb_clear(s);
    fmprb_clear(t);
    fmprb_clear(c);
}

void
fmprb_poly_evaluate2_rectangular(fmprb_t r, fmprb_t s, const fmprb_poly_t f, const fmprb_t a, long prec)
{
    _fmprb_poly_evaluate2_rectangular(r, s, f->coeffs, f->length, a, prec);
}

