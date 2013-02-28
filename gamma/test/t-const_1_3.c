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

#include "gamma.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("const_1_3....");
    fflush(stdout);
    flint_randinit(state);

    for (iter = 0; iter < 250; iter++)
    {
        fmprb_t r, s;
        fmpq_t q;
        long accuracy, prec;

        prec = 2 + n_randint(state, 1 << n_randint(state, 15));

        fmprb_init(r);
        fmprb_init(s);
        fmpq_init(q);

        fmpq_set_si(q, 1, 3);
        gamma_const_1_3(r, prec);
        gamma_series_fmpq_hypgeom(s, q, 1, prec);

        if (!fmprb_overlaps(r, s))
        {
            printf("FAIL: containment\n\n");
            printf("prec = %ld\n", prec);
            printf("r = "); fmprb_printd(r, prec / 3.33); printf("\n\n");
            abort();
        }

        accuracy = fmprb_rel_accuracy_bits(r);

        if (accuracy < prec - 4)
        {
            printf("FAIL: poor accuracy\n\n");
            printf("prec = %ld\n", prec);
            printf("r = "); fmprb_printd(r, prec / 3.33); printf("\n\n");
            abort();
        }

        fmprb_clear(r);
        fmprb_clear(s);
        fmpq_clear(q);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

