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

#include "fmpz_extras.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("lshift_mpn....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000; iter++)
    {
        fmpz_t a, b, c;
        ulong e;
        mp_limb_t atmp;
        mp_ptr aptr;
        mp_size_t an;
        int asgnbit;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);

        fmpz_randtest_not_zero(a, state, 1 + n_randint(state, 2000));
        fmpz_randtest(b, state, 1 + n_randint(state, 2000));
        fmpz_randtest(c, state, 1 + n_randint(state, 2000));

        e = n_randint(state, 1000);

        FMPZ_GET_MPN_READONLY(asgnbit, an, aptr, atmp, *a)
        fmpz_lshift_mpn(b, aptr, an, asgnbit, e);

        fmpz_mul_2exp(c, a, e);

        if (!fmpz_equal(b, c))
        {
            printf("FAIL\n");
            fmpz_print(a); printf("\n\n");
            fmpz_print(b); printf("\n\n");
            fmpz_print(c); printf("\n\n");
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

