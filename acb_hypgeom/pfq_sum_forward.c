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

    Copyright (C) 2014 Fredrik Johansson

******************************************************************************/

#include "acb_hypgeom.h"

void
acb_hypgeom_pfq_sum_forward(acb_t s, acb_t t,
    acb_srcptr a, long p, acb_srcptr b, long q, const acb_t z, long n, long prec)
{
    acb_t u, v;
    long k, i;

    acb_init(u);
    acb_init(v);

    acb_zero(s);
    acb_one(t);

    for (k = 0; k < n && !acb_is_zero(t); k++)
    {
        acb_add(s, s, t, prec);

        if (p > 0)
        {
            acb_add_ui(u, a, k, prec);

            for (i = 1; i < p; i++)
            {
                acb_add_ui(v, a + i, k, prec);
                acb_mul(u, u, v, prec);
            }

            acb_mul(t, t, u, prec);
        }

        if (q > 0)
        {
            acb_add_ui(u, b, k, prec);

            for (i = 1; i < q; i++)
            {
                acb_add_ui(v, b + i, k, prec);
                acb_mul(u, u, v, prec);
            }

            acb_div(t, t, u, prec);
        }

        acb_mul(t, t, z, prec);
    }

    acb_clear(u);
    acb_clear(v);
}

