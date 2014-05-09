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

#include "arb.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("get_unique_fmpz....");
    fflush(stdout);
    flint_randinit(state);

    for (iter = 0; iter < 100000; iter++)
    {
        arb_t x, z;
        fmpz_t y, a, b, exp;
        int unique, unique2;

        arb_init(x);
        arb_init(z);
        fmpz_init(y);
        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(exp);

        arb_randtest(x, state, 2000, 10);

        unique = arb_get_unique_fmpz(y, x);

        arb_get_interval_fmpz_2exp(a, b, exp, x);
        if (fmpz_sgn(exp) >= 0)
        {
            fmpz_mul_2exp(a, a, fmpz_get_si(exp));
            fmpz_mul_2exp(b, b, fmpz_get_si(exp));
        }
        else
        {
            fmpz_cdiv_q_2exp(a, a, -fmpz_get_si(exp));
            fmpz_fdiv_q_2exp(b, b, -fmpz_get_si(exp));
        }
        unique2 = fmpz_equal(a, b);

        if ((unique != unique2) || (unique && !fmpz_equal(y, a)))
        {
            printf("FAIL:\n\n");
            printf("x = "); arb_print(x); printf("\n\n");
            printf("unique = %d, unique2 = %d\n\n", unique, unique2);
            printf("y = "); fmpz_print(y); printf("\n\n");
            printf("a = "); fmpz_print(a); printf("\n\n");
            printf("b = "); fmpz_print(b); printf("\n\n");
            printf(" exp = "); fmpz_print(exp); printf("\n\n");
            abort();
        }

        if (unique)
        {
            arb_set_fmpz(z, y);
            arb_set_round(z, z, 2 + n_randint(state, 1000));

            if (!arb_overlaps(x, z))
            {
                printf("FAIL (overlap):\n\n");
                printf("x = "); arb_print(x); printf("\n\n");
                printf("y = "); fmpz_print(y); printf("\n\n");
                printf("z = "); arb_print(z); printf("\n\n");
                abort();
            }

            fmpz_add_ui(b, y, 1);
            if (arb_contains_fmpz(x, b))
            {
                printf("FAIL (contains a + 1):\n\n");
                printf("x = "); arb_print(x); printf("\n\n");
                printf("y = "); fmpz_print(y); printf("\n\n");
                abort();
            }

            fmpz_sub_ui(b, y, 1);
            if (arb_contains_fmpz(x, b))
            {
                printf("FAIL (contains a - 1):\n\n");
                printf("x = "); arb_print(x); printf("\n\n");
                printf("y = "); fmpz_print(y); printf("\n\n");
                abort();
            }
        }

        arb_clear(x);
        arb_clear(z);
        fmpz_clear(y);
        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(exp);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
