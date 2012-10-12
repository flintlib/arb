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
fmpz_holonomic_shift_fmpq(fmpz_holonomic_t res, const fmpz_holonomic_t op, const fmpq_t s)
{
    long i, j, d;

    if (fmpz_is_one(fmpq_denref(s)))
    {
        fmpz_holonomic_shift_fmpz(res, op, fmpq_numref(s));
    }
    else
    {
        fmpz_t Qpow;
        fmpz_poly_t spoly;
        fmpz C[2];

        if (res != op)
            fmpz_holonomic_set(res, op);

        d = fmpz_holonomic_degree(res);

        C[0] = *fmpq_numref(s);
        C[1] = *fmpq_denref(s);

        spoly->coeffs = C;
        spoly->length = 2;
        spoly->alloc = 2;

        fmpz_init(Qpow);
        fmpz_one(Qpow);

        for (j = d - 1; j >= 0; j--)
        {
            fmpz_mul(Qpow, Qpow, fmpq_denref(s));

            for (i = 0; i < res->length; i++)
            {
                fmpz_poly_struct * poly = res->coeffs + i;
                if (j < poly->length)
                    fmpz_mul(poly->coeffs + j, poly->coeffs + j, Qpow);
            }
        }

        for (i = 0; i < res->length; i++)
            fmpz_poly_compose(res->coeffs + i, res->coeffs + i, spoly);

        fmpz_clear(Qpow);
    }
}

