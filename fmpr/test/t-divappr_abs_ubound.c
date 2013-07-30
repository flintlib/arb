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

#include "fmpr.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("divappr_abs_ubound....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        fmpr_t a, b, c, d;
        long prec;

        fmpr_init(a);
        fmpr_init(b);
        fmpr_init(c);
        fmpr_init(d);

        fmpr_randtest_special(a, state, 2 + n_randint(state, 200), 100);
        fmpr_randtest_special(b, state, 2 + n_randint(state, 200), 100);
        fmpr_randtest_special(c, state, 2 + n_randint(state, 200), 100);
        fmpr_randtest_special(d, state, 2 + n_randint(state, 200), 100);
        prec = 2 + n_randint(state, 200);

        fmpr_div(c, a, b, prec, FMPR_RND_UP);
        fmpr_abs(c, c);

        fmpr_divappr_abs_ubound(d, a, b, prec);

        if (fmpr_cmp(c, d) > 0)
        {
            printf("FAIL:\n");
            fmpr_printd(a, prec / 3.32); printf("\n");
            fmpr_printd(b, prec / 3.32); printf("\n");
            fmpr_printd(c, prec / 3.32); printf("\n");
            fmpr_printd(d, prec / 3.32); printf("\n");
            abort();
        }

        fmpr_clear(a);
        fmpr_clear(b);
        fmpr_clear(c);
        fmpr_clear(d);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

