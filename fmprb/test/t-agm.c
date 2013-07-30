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

#include "fmprb.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("agm....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000; iter++)
    {
        fmprb_t a, b, c;
        fmpq_t q, r;
        mpfr_t t, u;
        long prec = 2 + n_randint(state, 200);

        fmprb_init(a);
        fmprb_init(b);
        fmprb_init(c);
        fmpq_init(q);
        fmpq_init(r);
        mpfr_init2(t, prec + 100);
        mpfr_init2(u, prec + 100);

        fmprb_randtest(a, state, 1 + n_randint(state, 200), 3);
        fmprb_randtest(b, state, 1 + n_randint(state, 200), 3);
        fmprb_randtest(c, state, 1 + n_randint(state, 200), 3);

        fmprb_agm(c, a, b, prec);

        if (fmprb_equal(a, b))
        {
            if (!fmprb_contains(c, a))
            {
                printf("FAIL: containment (identity)\n\n");
                printf("a = "); fmprb_print(a); printf("\n\n");
                printf("b = "); fmprb_print(b); printf("\n\n");
                printf("c = "); fmprb_print(c); printf("\n\n");
                abort();
            }
        }
        else
        {
            fmprb_get_rand_fmpq(q, state, a, 1 + n_randint(state, 200));
            fmprb_get_rand_fmpq(r, state, b, 1 + n_randint(state, 200));
            fmpq_get_mpfr(t, q, MPFR_RNDN);
            fmpq_get_mpfr(u, r, MPFR_RNDN);
            mpfr_agm(t, t, u, MPFR_RNDN);

            if (!fmprb_contains_mpfr(c, t))
            {
                printf("FAIL: containment\n\n");
                printf("a = "); fmprb_print(a); printf("\n\n");
                printf("b = "); fmprb_print(b); printf("\n\n");
                printf("c = "); fmprb_print(c); printf("\n\n");
                abort();
            }
        }

        fmprb_clear(a);
        fmprb_clear(b);
        fmprb_clear(c);
        fmpq_clear(q);
        fmpq_clear(r);
        mpfr_clear(t);
        mpfr_clear(u);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

