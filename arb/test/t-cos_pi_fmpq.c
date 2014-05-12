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

#include "arb.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("cos_pi_fmpq....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        arb_t c1, c2;
        fmpq_t x;
        long prec;

        prec = 2 + n_randint(state, 5000);

        arb_init(c1);
        arb_init(c2);
        fmpq_init(x);

        fmpq_randtest(x, state, 1 + n_randint(state, 200));

        arb_cos_pi_fmpq(c1, x, prec);

        arb_const_pi(c2, prec);
        arb_mul_fmpz(c2, c2, fmpq_numref(x), prec);
        arb_div_fmpz(c2, c2, fmpq_denref(x), prec);
        arb_cos(c2, c2, prec);

        if (!arb_overlaps(c1, c2))
        {
            printf("FAIL: overlap\n\n");
            printf("x = "); fmpq_print(x); printf("\n\n");
            printf("c1 = "); arb_printd(c1, 15); printf("\n\n");
            printf("c2 = "); arb_printd(c2, 15); printf("\n\n");
            abort();
        }

        arb_clear(c1);
        arb_clear(c2);
        fmpq_clear(x);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

