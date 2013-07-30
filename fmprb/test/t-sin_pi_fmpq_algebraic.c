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

    printf("sin_pi_fmpq_algebraic....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        fmprb_t s1, s2;
        ulong p, q, g;
        long prec;

        prec = 2 + n_randint(state, 5000);
        q = 1 + n_randint(state, 500);
        p = n_randint(state, q / 2 + 1);

        g = n_gcd(q, p);
        q /= g;
        p /= g;

        fmprb_init(s1);
        fmprb_init(s2);

        _fmprb_sin_pi_fmpq_algebraic(s1, p, q, prec);

        fmprb_const_pi(s2, prec);
        fmprb_mul_ui(s2, s2, p, prec);
        fmprb_div_ui(s2, s2, q, prec);
        fmprb_sin(s2, s2, prec);

        if (!fmprb_overlaps(s1, s2))
        {
            printf("FAIL: overlap\n\n");
            printf("p/q = %lu/%lu", p, q); printf("\n\n");
            printf("s1 = "); fmprb_printd(s1, 15); printf("\n\n");
            printf("s2 = "); fmprb_printd(s2, 15); printf("\n\n");
            abort();
        }

        if (fmprb_rel_accuracy_bits(s1) < prec - 2)
        {
            printf("FAIL: accuracy\n\n");
            printf("p/q = %lu/%lu", p, q); printf("\n\n");
            printf("prec=%ld eff=%ld\n", prec, fmprb_rel_accuracy_bits(s1));
            printf("s1 = "); fmprb_printd(s1, 15); printf("\n\n");
            printf("s2 = "); fmprb_printd(s2, 15); printf("\n\n");
            abort();
        }

        fmprb_clear(s1);
        fmprb_clear(s2);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

