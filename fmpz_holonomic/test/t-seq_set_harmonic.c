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
#include "arith.h"

int main()
{
    printf("seq_set_harmonic....");
    fflush(stdout);

    {
        long i;
        fmpq_t t, u;
        fmpz_holonomic_t op;
        fmpq initial[2];

        fmpz_holonomic_init(op);
        fmpq_init(t);
        fmpq_init(u);

        fmpz_holonomic_seq_set_harmonic(op);

        *fmpq_numref(&initial[0]) = 0;
        *fmpq_denref(&initial[0]) = 1;

        *fmpq_numref(&initial[1]) = 1;
        *fmpq_denref(&initial[1]) = 1;

        for (i = 0; i < 20; i++)
        {
            fmpz_holonomic_get_nth_fmpq(t, op, initial, 0, i);
            arith_harmonic_number(u, i);
            if (!fmpq_equal(t, u))
            {
                printf("FAIL\n");
                printf("i = %ld, t = ", i); fmpq_print(t);
                printf("   u = "); fmpq_print(u);
                printf("\n");
                abort();
            }
        }

        fmpq_clear(t);
        fmpq_clear(u);
        fmpz_holonomic_clear(op);

    }

    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

