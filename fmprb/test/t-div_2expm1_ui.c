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

#include "fmprb.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("div_2expm1_ui....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000; iter++)
    {
        fmprb_t a, b, c;
        ulong n;
        long prec, acc1, acc2;
        fmpz_t t;

        fmprb_init(a);
        fmprb_init(b);
        fmprb_init(c);
        fmpz_init(t);

        prec = 2 + n_randint(state, 10000);
        fmprb_randtest(a, state, 1 + n_randint(state, 10000), 10);
        n = 1 + (n_randtest(state) % (10 * prec));

        fmprb_div_2expm1_ui(b, a, n, prec);

        fmprb_one(c);
        fmpz_set_ui(fmpr_expref(fmprb_midref(c)), n);
        fmprb_sub_ui(c, c, 1, prec);
        fmprb_div(c, a, c, prec);

        acc1 = fmprb_rel_accuracy_bits(a);
        acc2 = fmprb_rel_accuracy_bits(b);

        if (!fmprb_overlaps(b, c))
        {
            printf("FAIL: containment\n\n");
            printf("a = "); fmprb_print(a); printf("\n\n");
            printf("b = "); fmprb_print(b); printf("\n\n");
            printf("c = "); fmprb_print(c); printf("\n\n");
            abort();
        }

        if ((acc2 < FLINT_MIN(prec, acc1) - 10) &&
                !(acc1 == -FMPR_PREC_EXACT && acc2 == -FMPR_PREC_EXACT))
        {
            printf("FAIL: poor accuracy\n\n");
            printf("prec=%ld, acc1=%ld, acc2=%ld\n\n", prec, acc1, acc2);
            printf("n = %lu\n\n", n);
            printf("a = "); fmprb_print(a); printf("\n\n");
            printf("b = "); fmprb_print(b); printf("\n\n");
            printf("c = "); fmprb_print(c); printf("\n\n");
            abort();
        }

        fmprb_clear(a);
        fmprb_clear(b);
        fmprb_clear(c);
        fmpz_clear(t);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

