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

    Copyright (C) 2015 Fredrik Johansson

******************************************************************************/

#include "acb_hypgeom.h"

void
acb_hypgeom_legendre_p(acb_t res, const acb_t n, const acb_t m,
    const acb_t z, int type, long prec)
{
    acb_t a, b, c, w;

    acb_init(a);
    acb_init(b);
    acb_init(c);
    acb_init(w);

    acb_neg(a, n);
    acb_add_ui(b, n, 1, prec);

    acb_sub_ui(c, m, 1, prec);
    acb_neg(c, c);

    acb_sub_ui(w, z, 1, prec);
    acb_neg(w, w);
    acb_mul_2exp_si(w, w, -1);

    acb_hypgeom_2f1(w, a, b, c, w, 1, prec);

    if (!acb_is_zero(m))
    {
        acb_add_ui(a, z, 1, prec);
        acb_sub_ui(b, z, 1, prec);

        if (type == 0)
        {
            acb_neg(b, b);
        }
        else if (type != 1)
        {
            printf("unsupported 'type' %d for legendre p\n", type);
            abort();
        }

        acb_mul_2exp_si(c, m, -1);
        acb_pow(a, a, c, prec);
        acb_neg(c, c);
        acb_pow(b, b, c, prec);
        acb_mul(w, w, a, prec);
        acb_mul(w, w, b, prec);
    }

    acb_set(res, w);

    acb_clear(a);
    acb_clear(b);
    acb_clear(c);
    acb_clear(w);
}

