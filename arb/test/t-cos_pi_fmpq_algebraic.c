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

    printf("cos_pi_fmpq_algebraic....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        arb_t c1, c2;
        ulong p, q, g;
        long prec;

        prec = 2 + n_randint(state, 5000);
        q = 1 + n_randint(state, 500);
        p = n_randint(state, q / 2 + 1);

        g = n_gcd(q, p);
        q /= g;
        p /= g;

        arb_init(c1);
        arb_init(c2);

        _arb_cos_pi_fmpq_algebraic(c1, p, q, prec);

        arb_const_pi(c2, prec);
        arb_mul_ui(c2, c2, p, prec);
        arb_div_ui(c2, c2, q, prec);
        arb_cos(c2, c2, prec);

        if (!arb_overlaps(c1, c2))
        {
            printf("FAIL: overlap\n\n");
            printf("p/q = %lu/%lu", p, q); printf("\n\n");
            printf("c1 = "); arb_printd(c1, 15); printf("\n\n");
            printf("c2 = "); arb_printd(c2, 15); printf("\n\n");
            abort();
        }

        if (arb_rel_accuracy_bits(c1) < prec - 2)
        {
            printf("FAIL: accuracy\n\n");
            printf("p/q = %lu/%lu", p, q); printf("\n\n");
            printf("prec=%ld eff=%ld\n", prec, arb_rel_accuracy_bits(c1));
            printf("c1 = "); arb_printd(c1, 15); printf("\n\n");
            printf("c2 = "); arb_printd(c2, 15); printf("\n\n");
            abort();
        }

        arb_clear(c1);
        arb_clear(c2);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

