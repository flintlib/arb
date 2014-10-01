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

#include "acb_modular.h"

void
psl2z_randtest(psl2z_t g, flint_rand_t state, long bits)
{
    bits = FLINT_MAX(bits, 1);

    fmpz_randtest(&g->a, state, bits);
    fmpz_randtest(&g->b, state, bits);

    if (fmpz_is_zero(&g->a) && fmpz_is_zero(&g->b))
    {
        psl2z_one(g);
    }
    else
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_xgcd(t, &g->d, &g->c, &g->a, &g->b);
        fmpz_divexact(&g->a, &g->a, t);
        fmpz_divexact(&g->b, &g->b, t);

        if (fmpz_sgn(&g->c) < 0)
            fmpz_neg(&g->c, &g->c);
        else
            fmpz_neg(&g->b, &g->b);

        if (fmpz_is_zero(&g->c) && fmpz_sgn(&g->d) < 0)
        {
            fmpz_neg(&g->a, &g->a);
            fmpz_neg(&g->b, &g->b);
            fmpz_neg(&g->c, &g->c);
            fmpz_neg(&g->d, &g->d);
        }

        fmpz_clear(t);
    }
}

