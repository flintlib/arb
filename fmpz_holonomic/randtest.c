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

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "fmpz_holonomic.h"

void
fmpz_holonomic_randtest(fmpz_holonomic_t op, flint_rand_t state, long r, long d, long b)
{
    long i;

    r = n_randint(state, r + 1);
    d = n_randint(state, d + 1);
    b = n_randint(state, b + 1);

    fmpz_holonomic_fit_length(op, r + 1);

    for (i = 0; i <= r; i++)
        fmpz_poly_randtest(op->coeffs + i, state, d + 1, b);

    _fmpz_holonomic_set_length(op, r + 1);
    fmpz_holonomic_normalise_leading(op);

    if (op->length == 0)
        fmpz_holonomic_one(op);
}

