/*=============================================================================

    This file is part of ARB.

    ARB is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    ARB is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with ARB; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "arb.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("exp....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000; iter++)
    {
        arb_t r, s;
        fmpq_t q;
        mpfr_t x, y;
        long wp;

        arb_init(r, 1 + n_randint(state, 600));
        arb_init(s, 1 + n_randint(state, 600));

        wp = FLINT_MAX(arb_prec(r), arb_prec(s)) + 100;

        mpfr_init2(x, wp);
        mpfr_init2(y, wp);

        fmpq_init(q);

        do {
            arb_randtest(r, state, 5);
        } while (FLINT_MAX(fmpz_bits(arb_midref(r)),
            fmpz_bits(arb_radref(r))) + *arb_expref(r) > 13);

        arb_get_rand_fmpq(q, state, r);
        fmpq_get_mpfr(x, q, MPFR_RNDN);

        arb_exp(s, r);
        mpfr_exp(y, x, MPFR_RNDN);

        if (!arb_contains_mpfr(s, y))
        {
            printf("FAIL: containment\n\n");
            printf("r = "); arb_debug(r); printf("\n\n");
            printf("s = "); arb_debug(s); printf("\n\n");
            abort();
        }

        arb_clear(r);
        arb_clear(s);

        fmpq_clear(q);

        mpfr_clear(x);
        mpfr_clear(y);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    mpfr_free_cache();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
