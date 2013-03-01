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

    printf("small_frac....");
    fflush(stdout);
    flint_randinit(state);

    for (iter = 0; iter < 1000; iter++)
    {
        fmprb_t r, s;
        fmpq_t q;
        long accuracy, prec;

        prec = 2 + n_randint(state, 1 << n_randint(state, 15));

        fmprb_init(r);
        fmprb_init(s);
        fmpq_init(q);

        switch (n_randint(state, 8))
        {
            case 0: fmpq_set_si(q, 1, 1); break;
            case 1: fmpq_set_si(q, 1, 2); break;
            case 2: fmpq_set_si(q, 1, 3); break;
            case 3: fmpq_set_si(q, 2, 3); break;
            case 4: fmpq_set_si(q, 1, 4); break;
            case 5: fmpq_set_si(q, 3, 4); break;
            case 6: fmpq_set_si(q, 1, 6); break;
            default: fmpq_set_si(q, 5, 6); break;
        }

        gamma_small_frac(r, *fmpq_numref(q), *fmpq_denref(q), prec);
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

