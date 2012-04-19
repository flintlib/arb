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

    printf("zeta_ui_euler_product....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        arb_t r;
        ulong n;
        mpfr_t s;
        long effective_prec;

        arb_init(r, 1 + n_randint(state, 1 << n_randint(state, 14)));
        mpfr_init2(s, arb_prec(r) + 100);
        arb_randtest(r, state, 10);

        /* don't take too small arguments (for fast convergence) */
        n = FLINT_MAX(6, (0.06 * arb_prec(r))) +
            2*n_randint(state, arb_prec(r));

        arb_zeta_inv_ui_euler_product(r, n);
        mpfr_zeta_ui(s, n, MPFR_RNDN);
        mpfr_ui_div(s, 1, s, MPFR_RNDN);

        if (!arb_contains_mpfr(r, s))
        {
            printf("FAIL: containment\n\n");
            printf("n = %lu\n\n", n);
            printf("r = "); arb_debug(r); printf("\n\n");
            printf("s = "); mpfr_printf("%.275Rf\n", s); printf("\n\n");
            abort();
        }

        effective_prec = fmpz_bits(arb_midref(r)) - fmpz_bits(arb_radref(r));

        if (!fmpz_is_zero(arb_radref(r)) && effective_prec < arb_prec(r) - 4)
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
