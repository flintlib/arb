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

    printf("const_zeta3_bsplit....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 250; iter++)
    {
        arb_t r;
        mpfr_t s;
        long effective_prec;

        arb_init(r, 1 + n_randint(state, 1 << n_randint(state, 17)));
        mpfr_init2(s, arb_prec(r) + 1000);
        arb_randtest(r, state, 10);

        arb_const_zeta3_bsplit(r);
        mpfr_zeta_ui(s, 3, MPFR_RNDN);

        if (!arb_contains_mpfr(r, s))
        {
            printf("FAIL: containment\n\n");
            printf("r = "); arb_debug(r); printf("\n\n");
            abort();
        }

        effective_prec = fmpz_bits(arb_midref(r)) - fmpz_bits(arb_radref(r));

        if (effective_prec < arb_prec(r) - 4)
        {
            printf("FAIL: poor accuracy\n\n");
            printf("r = "); arb_debug(r); printf("\n\n");
            abort();
        }

        arb_clear(r);
        mpfr_clear(s);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    mpfr_free_cache();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
