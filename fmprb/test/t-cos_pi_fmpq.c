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

    printf("cos_pi_fmpq....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        fmprb_t c1, c2;
        fmpq_t x;
        long prec;

        prec = 2 + n_randint(state, 1000);

        fmprb_init(c1);
        fmprb_init(c2);
        fmpq_init(x);

        fmpq_randtest(x, state, 1 + n_randint(state, 200));

        fmprb_cos_pi_fmpq(c1, x, prec);

        fmprb_const_pi(c2, prec);
        fmprb_mul_fmpz(c2, c2, fmpq_numref(x), prec);
        fmprb_div_fmpz(c2, c2, fmpq_denref(x), prec);
        fmprb_cos(c2, c2, prec);

        if (!fmprb_overlaps(c1, c2))
        {
            printf("FAIL: overlap\n\n");
            printf("x = "); fmpq_print(x); printf("\n\n");
            printf("c1 = "); fmprb_printd(c1, 15); printf("\n\n");
            printf("c2 = "); fmprb_printd(c2, 15); printf("\n\n");
            abort();
        }

        fmprb_clear(c1);
        fmprb_clear(c2);
        fmpq_clear(x);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

