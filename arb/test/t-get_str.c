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

    Copyright (C) 2015 Fredrik Johansson

******************************************************************************/

#include <string.h>
#include "arb.h"

int main()
{
    flint_rand_t state;
    long iter;

    printf("get_str....");
    fflush(stdout);
    flint_randinit(state);

    /* just test no crashing... */
    for (iter = 0; iter < 10000; iter++)
    {
        arb_t x;
        char * s;
        long n;

        arb_init(x);

        arb_randtest_special(x, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));

        n = 1 + n_randint(state, 300);

        s = arb_get_str(x, n, (n_randint(state, 2) * ARB_STR_MORE)
                            |  (n_randint(state, 2) * ARB_STR_NO_RADIUS)
                            | (ARB_STR_CONDENSE * n_randint(state, 50)));

        flint_free(s);
        arb_clear(x);
    }

    for (iter = 0; iter < 100000; iter++)
    {
        arb_t x, y;
        char * s;
        long n, prec;
        int conversion_error;

        arb_init(x);
        arb_init(y);

        arb_randtest_special(x, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        arb_randtest_special(y, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));

        n = 1 + n_randint(state, 300);
        prec = 2 + n_randint(state, 1000);

        s = arb_get_str(x, n, n_randint(state, 2) * ARB_STR_MORE);
        conversion_error = arb_set_str(y, s, prec);

        if (conversion_error || !arb_contains(y, x))
        {
            printf("FAIL (roundtrip)  iter = %ld\n", iter);
            printf("x = "); arb_printd(x, 50); printf("\n\n");
            printf("s = %s", s); printf("\n\n");
            printf("y = "); arb_printd(y, 50); printf("\n\n");
            abort();
        }

        flint_free(s);
        arb_clear(x);
        arb_clear(y);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

