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

    printf("gamma_fmpq....");
    fflush(stdout);
    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        arb_t r, s;
        fmpq_t q;
        long accuracy, prec, pp, qq;

        prec = 2 + n_randint(state, 1 << n_randint(state, 12));
        prec += 20;

        arb_init(r);
        arb_init(s);
        fmpq_init(q);

        pp = -100 + n_randint(state, 10000);
        qq = 1 + n_randint(state, 20);
        fmpq_set_si(q, pp, qq);

        arb_gamma_fmpq(r, q, prec);

        arb_set_fmpq(s, q, prec);
        arb_gamma(s, s, prec);

        if (!arb_overlaps(r, s))
        {
            printf("FAIL: containment\n\n");
            printf("prec = %ld\n", prec);
            printf("q = "); fmpq_print(q); printf("\n\n");
            printf("r = "); arb_printd(r, prec / 3.33); printf("\n\n");
            printf("s = "); arb_printd(s, prec / 3.33); printf("\n\n");
            abort();
        }

        if (!(fmpz_is_one(fmpq_denref(q)) && fmpz_sgn(fmpq_numref(q)) <= 0)
            && FLINT_ABS(pp / qq) < 10)
        {
            accuracy = arb_rel_accuracy_bits(r);

            if (accuracy < prec - 6)
            {
                printf("FAIL: poor accuracy\n\n");
                printf("prec = %ld\n", prec);
                printf("q = "); fmpq_print(q); printf("\n\n");
                printf("r = "); arb_printd(r, prec / 3.33); printf("\n\n");
                abort();
            }
        }

        arb_clear(r);
        arb_clear(s);
        fmpq_clear(q);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

