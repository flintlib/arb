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

    printf("stirling_eval_fmpcb_series....");
    fflush(stdout);

    flint_randinit(state);

    /* check self-consistency */
    for (iter = 0; iter < 1000; iter++)
    {
        fmpcb_ptr u, v;
        fmpcb_t x;
        long i, n, len1, len2, prec1, prec2;

        len1 = 1 + n_randint(state, 40);
        len2 = 1 + n_randint(state, 40);

        prec1 = 2 + n_randint(state, 1000);
        prec2 = 2 + n_randint(state, 1000);

        n = 1 + n_randint(state, 40);

        u = _fmpcb_vec_init(len1);
        v = _fmpcb_vec_init(len2);
        fmpcb_init(x);

        fmpcb_randtest(x, state, 2 + n_randint(state, 1000), 20);

        gamma_stirling_eval_fmpcb_series(u, x, n, len1, prec1);
        gamma_stirling_eval_fmpcb_series(v, x, n, len2, prec2);

        for (i = 0; i < FLINT_MIN(len1, len2); i++)
        {
            if (!fmpcb_overlaps(u + i, v + i))
            {
                printf("FAIL: overlap\n\n");
                printf("n = %ld, len1 = %ld, len2 = %ld, i = %ld\n\n", n, len1, len2, i);
                printf("x = "); fmpcb_printd(x, prec1 / 3.33); printf("\n\n");
                printf("u = "); fmpcb_printd(u + i, prec1 / 3.33); printf("\n\n");
                printf("v = "); fmpcb_printd(v + i, prec2 / 3.33); printf("\n\n");
                abort();
            }
        }

        _fmpcb_vec_clear(u, len1);
        _fmpcb_vec_clear(v, len2);
        fmpcb_clear(x);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

