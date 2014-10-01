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

int
psl2z_is_correct(const psl2z_t g)
{
    int res;
    fmpz_t t;

    if (fmpz_sgn(&g->c) < 0)
        return 0;

    if (fmpz_is_zero(&g->c) && fmpz_sgn(&g->d) <= 0)
        return 0;

    fmpz_init(t);
    fmpz_mul(t, &g->a, &g->d);
    fmpz_submul(t, &g->b, &g->c);
    res = fmpz_is_one(t);
    fmpz_clear(t);
    return res;
}

