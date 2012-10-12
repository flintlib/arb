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
fmpz_holonomic_get_nth_fmpq(fmpq_t res, const fmpz_holonomic_t op, const fmpq * initial, long n0, long n)
{
    long r = fmpz_holonomic_order(op);

    if (r == 0)
    {
        fmpq_zero(res);
        return;
    }
    else if (n < n0)
    {
        printf("not implemented\n");
        abort();
    }
    else if (n - n0 < r)
    {
        fmpq_set(res, initial + n - n0);
        return;
    }
    else
    {
        fmpz_mat_t M;
        long i;
        fmpz_t Q;
        fmpz_mat_init(M, r, r);
        fmpz_init(Q);

        fmpz_holonomic_forward_fmpz_mat(M, Q, op, n0, n - n0 - r + 1);

        {
            fmpz_t g, t;
            fmpz_init(g);
            fmpz_init(t);

            fmpz_one(g);
            for (i = 0; i < r; i++)
                fmpz_lcm(g, g, fmpq_denref(initial + i));

            fmpz_divexact(t, g, fmpq_denref(initial + 0));
            fmpz_mul(t, t, fmpq_numref(initial + 0));
            fmpz_mul(fmpz_mat_entry(M, r - 1, 0), fmpz_mat_entry(M, r - 1, 0), t);

            for (i = 1; i < r; i++)
            {
                fmpz_divexact(t, g, fmpq_denref(initial + i));
                fmpz_mul(t, t, fmpq_numref(initial + i));
                fmpz_addmul(fmpz_mat_entry(M, r - 1, 0), fmpz_mat_entry(M, r - 1, i), t);
            }

            fmpz_set(fmpq_numref(res), fmpz_mat_entry(M, r - 1, 0));
            fmpz_mul(fmpq_denref(res), Q, g);
            fmpq_canonicalise(res);

            fmpz_clear(g);
            fmpz_clear(t);
        }

        fmpz_mat_clear(M);
        fmpz_clear(Q);
    }
}

