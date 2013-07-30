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

    Copyright (C) 2013 Fredrik Johansson

******************************************************************************/

#include "fmpz_extras.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("sub_si_inline....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000; iter++)
    {
        fmpz_t a, c, d;
        long b;

        fmpz_init(a);
        fmpz_init(c);
        fmpz_init(d);

        fmpz_randtest(a, state, 1 + n_randint(state, 200));
        fmpz_randtest(c, state, 1 + n_randint(state, 200));
        fmpz_randtest(d, state, 1 + n_randint(state, 200));
        b = n_randtest(state);

        fmpz_sub_si(c, a, b);
        fmpz_sub_si_inline(d, a, b);

        if (!fmpz_equal(c, d))
        {
            printf("FAIL\n");
            fmpz_print(a); printf("\n\n");
            printf("%ld", b); printf("\n\n");
            fmpz_print(c); printf("\n\n");
            fmpz_print(d); printf("\n\n");
            abort();
        }

        fmpz_sub_si_inline(a, a, b);
        if (!fmpz_equal(c, a))
        {
            printf("FAIL (aliasing)\n");
            fmpz_print(a); printf("\n\n");
            printf("%ld", b); printf("\n\n");
            fmpz_print(c); printf("\n\n");
            fmpz_print(d); printf("\n\n");
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(c);
        fmpz_clear(d);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

