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
psl2z_inv(psl2z_t h, const psl2z_t g)
{
    if (h != g)
        psl2z_set(h, g);

    if (fmpz_is_zero(&h->c) && fmpz_sgn(&h->a) > 0)
    {
        fmpz_neg(&h->b, &h->b);
        fmpz_swap(&h->d, &h->a);
    }
    else
    {
        fmpz_swap(&h->a, &h->d);
        fmpz_neg(&h->a, &h->a);
        fmpz_neg(&h->d, &h->d);
    }
}

