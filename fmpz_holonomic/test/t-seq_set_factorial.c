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

int main()
{
    printf("seq_set_factorial....");
    fflush(stdout);

    {
        long i;
        fmpz_t t, u;
        fmpz_holonomic_t op;
        fmpz initial[1];

        fmpz_holonomic_init(op);
        fmpz_init(t);
        fmpz_init(u);

        fmpz_holonomic_seq_set_factorial(op);

        initial[0] = 1;

        for (i = 0; i < 20; i++)
        {
            fmpz_holonomic_get_nth_fmpz(t, op, initial, 0, i);
            fmpz_fac_ui(u, i);
            if (!fmpz_equal(t, u))
            {
                printf("FAIL\n");
                printf("i = %ld, t = ", i); fmpz_print(t);
                printf("   u = "); fmpz_print(u);
                printf("\n");
                abort();
            }
        }

        fmpz_clear(t);
        fmpz_clear(u);
        fmpz_holonomic_clear(op);

    }

    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

