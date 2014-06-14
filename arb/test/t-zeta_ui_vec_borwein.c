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

    printf("zeta_ui_vec_borwein....");
    fflush(stdout);
    flint_randinit(state);

    for (iter = 0; iter < 100; iter++)
    {
        arb_ptr r;
        ulong n;
        long i, num, step;
        mpfr_t s;
        long prec, accuracy;

        prec = 2 + n_randint(state, 1 << n_randint(state, 13));
        num = 1 + n_randint(state, 20);
        step = 1 + n_randint(state, 5);

        r = _arb_vec_init(num);
        mpfr_init2(s, prec + 100);

        do { n = n_randint(state, 1 << n_randint(state, 10)); } while (n < 2);

        arb_zeta_ui_vec_borwein(r, n, num, step, prec);

        for (i = 0; i < num; i++)
        {
            mpfr_zeta_ui(s, n + i * step, MPFR_RNDN);

            if (!arb_contains_mpfr(r + i, s))
            {
                printf("FAIL: containment\n\n");
                printf("n = %lu\n\n", n + i * step);
                printf("r = "); arb_printd(r + i, prec / 3.33); printf("\n\n");
                printf("s = "); mpfr_printf("%.275Rf\n", s); printf("\n\n");
                abort();
            }

            accuracy = arb_rel_accuracy_bits(r + i);

            if (accuracy < prec - 4)
            {
                printf("FAIL: accuracy = %ld, prec = %ld\n\n", accuracy, prec);
                printf("n = %lu\n\n", n + i * step);
                printf("r = "); arb_printd(r + i, prec / 3.33); printf("\n\n");
                abort();
            }
        }

        _arb_vec_clear(r, num);
        mpfr_clear(s);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
