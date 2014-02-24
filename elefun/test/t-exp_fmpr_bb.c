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

#include "elefun.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("exp_fmpr_bb....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        fmprb_t x, y, z;
        long prec = 2 + n_randint(state, 5000);

        fmprb_init(x);
        fmprb_init(y);
        fmprb_init(z);

        fmprb_randtest(x, state, 1 + n_randint(state, 5000), 3);
        fmpr_zero(fmprb_radref(x));

        if (n_randint(state, 2))
        {
#if FLINT_BITS == 64
            fmprb_mul_2exp_si(x, x, 1 + n_randint(state, 40));
#else
            fmprb_mul_2exp_si(x, x, 1 + n_randint(state, 20));
#endif
        }
        else
        {
            fmprb_mul_2exp_si(x, x, -n_randint(state, 1.5 * prec));
        }

        elefun_exp_via_mpfr(y, x, prec + 100);
        elefun_exp_fmpr_bb(z, fmprb_midref(x), prec, 0);

        if (!fmprb_contains(z, y))
        {
            printf("FAIL: containment\n\n");
            printf("prec = %ld\n\n", prec);
            printf("x = "); fmprb_print(x); printf("\n\n");
            printf("y = "); fmprb_print(y); printf("\n\n");
            printf("z = "); fmprb_print(z); printf("\n\n");
            abort();
        }

        if (fmprb_rel_accuracy_bits(z) < prec - 2)
        {
            printf("FAIL: poor accuracy\n\n");
            printf("prec = %ld,  acc = %ld\n\n", prec, fmprb_rel_accuracy_bits(z));
            printf("x = "); fmprb_print(x); printf("\n\n");
            printf("y = "); fmprb_print(y); printf("\n\n");
            printf("z = "); fmprb_print(z); printf("\n\n");
            abort();
        }

        elefun_exp_fmpr_bb(x, fmprb_midref(x), prec, 0);

        if (!fmprb_overlaps(x, z))
        {
            printf("FAIL: aliasing\n\n");
            abort();
        }

        fmprb_clear(x);
        fmprb_clear(y);
        fmprb_clear(z);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

