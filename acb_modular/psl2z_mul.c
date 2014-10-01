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
psl2z_mul(psl2z_t h, const psl2z_t f, const psl2z_t g)
{
    if (h == f || h == g)
    {
        psl2z_t t;
        psl2z_init(t);
        psl2z_mul(t, f, g);
        psl2z_swap(t, h);
        psl2z_clear(t);
        return;
    }

    fmpz_mul(&h->a, &f->a, &g->a);
    fmpz_addmul(&h->a, &f->b, &g->c);

    fmpz_mul(&h->b, &f->a, &g->b);
    fmpz_addmul(&h->b, &f->b, &g->d);

    fmpz_mul(&h->c, &f->c, &g->a);
    fmpz_addmul(&h->c, &f->d, &g->c);

    fmpz_mul(&h->d, &f->c, &g->b);
    fmpz_addmul(&h->d, &f->d, &g->d);

    if (fmpz_sgn(&h->c) < 0 || (fmpz_is_zero(&h->c) && fmpz_sgn(&h->d) < 0))
    {
        fmpz_neg(&h->a, &h->a);
        fmpz_neg(&h->b, &h->b);
        fmpz_neg(&h->c, &h->c);
        fmpz_neg(&h->d, &h->d);
    }
}

